"""
Tetracorder operations reimplemented in numpy / scipy.

The feature-fitting functions reproduce the native Tetracorder 6.00 /
specpr code paths:

    linear_feature_continuum   bdmset (reference) / bandmp (target)
    curved_feature_continuum   bandmpcv
    characterise_feature       bandmp / bandmpcv least squares, plus the
                               fit checks at the top of tp1mat
    fuzzy_greater / fuzzy_less tp1mat fuzzy logic

Native behaviour that these functions follow:

  * Window wavelengths are turned into channel ranges with specpr wtochbin
    on the full wavelength array.
  * tetracorder.r sets imgflg = 0 (since Tetracorder 3.1, 1995), so the
    continuum-removed spectrum and the least-squares sums run over BOTH
    continuum windows and the band (cl1..cr2), not only cl2+1..cr1-1.
  * The band minimum/maximum and the feature weight area use only the
    channels between the continuum windows.
  * A channel is skipped wherever the (continuum-removed) reference is
    deleted. Pass valid_bands & isfinite(reference) for the target so the
    target skips the same channels.
  * REAL*4 arithmetic, with the least-squares sums in REAL*8.
"""

from dataclasses import dataclass, field
from typing import Literal
from copy import deepcopy
from pathlib import Path
import sqlite3

import numpy as np
from scipy.interpolate import CubicSpline

from .config import VARIABLE_PRESETS


_TINY = np.float32(0.1e-20)       # specpr "effectively zero"
_CONT_MIN = 0.1e-6                # tp1mat default ct/lct/rct limits
_CONT_MAX = 0.1e20


def prepare_rules(rules, mode="default"):
    """
    Expects rules from the parsed json rules json
    mode = default, default_bak23, emit_c, MMM_09c, MMM255t,
    HYB2ryug defined in Tetracorder 6.00a VARIABLES/cmd.lib.setup.variables-*
    Custom modes can be established by creating a new preset dictionary
    """
    values = VARIABLE_PRESETS[mode]

    def replace(value):

        if isinstance(value, dict):
            return {key: replace(item) for key, item in value.items()}

        if isinstance(value, list):
            # Parser stored symbolic values as one-item lists.
            # Collapse the list and substitute the symbol directly.
            if (len(value) == 1 and isinstance(value[0], str) and value[0] in values):
                return deepcopy(values[value[0]])

            return [replace(item) for item in value]

        if isinstance(value, str) and value in values:
            return deepcopy(values[value])

        return value

    return replace(rules)


#========= load the prepared db of references ================================

def load_references(db_file: str | Path) -> dict:
    """
    Access the prepared SQLite database of reference spectra used in the Tetracorder
    algorithm.
    Returns a dictionary keyed by (library, record number) for use by the evaluator
    """
    references = {}

    connection = sqlite3.connect(db_file)
    connection.row_factory = sqlite3.Row

    rows = connection.execute(
        """
        SELECT
            Samples.SampleID,
            Samples.Name,
            Samples.Library,
            Samples.NativeRecord,
            Samples.ConvolvedRecord,
            Samples.OriginalTitle,
            Samples.InstrumentCode,
            Samples.PurityCode,
            Samples.MeasurementCode,
            Samples.SourceFilename,
            Spectra.XData,
            Spectra.YData,
            Spectra.FData
        FROM Samples
        JOIN Spectra USING (SampleID)
        WHERE Samples.ConvolvedRecord IS NOT NULL
        """
    )

    for row in rows:
        key = (row["Library"], row["ConvolvedRecord"])

        references[key] = {
            "sample_id": row["SampleID"],
            "name": row["Name"],
            "library": row["Library"],
            "native_record": row["NativeRecord"],
            "convolved_record": row["ConvolvedRecord"],
            "original_title": row["OriginalTitle"],
            "instrument_code": row["InstrumentCode"],
            "purity_code": row["PurityCode"],
            "measurement_code": row["MeasurementCode"],
            "source_filename": row["SourceFilename"],
            "wavelengths": np.frombuffer(row["XData"], dtype=np.float32).copy(),
            "reflectance": np.frombuffer(row["YData"], dtype=np.float32).copy(),
            "fwhm": np.frombuffer(row["FData"], dtype=np.float32).copy(),
        }

    connection.close()

    return references


#========= Helpers for the rule evaluation logic ==============================

def fuzzy_greater(value, thresholds):
    """
    tp1mat fuzzy logic for ">" tests and material constraints:

        x < z1 -> 0;  elif x < z2 -> (x - z1) / (z2 - z1);  else 1

    Limits are used in the order written, so a "backwards" pair
    (z1 > z2) is a hard cut at z1, as in the native code.
    NaN compares False everywhere and so gives 1, as in Fortran.
    The input dtype is kept (float32 in -> REAL*4 arithmetic).
    """
    z1, z2 = (float(t) for t in thresholds)
    x = np.asarray(np.ma.filled(np.ma.asarray(value), np.nan))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(x < z1, 0.0,
                        np.where(x < z2, (x - z1) / (z2 - z1), 1.0))


def fuzzy_less(value, thresholds):
    """
    tp1mat fuzzy logic for "<" tests (rcbblc<, lcbbrc<), limits in the
    order written:

        x > z1 -> 0;  elif x > z2 -> (x - z1) / (z2 - z1);  else 1

    e.g. "rcbblc< 0.8 0.9" is a hard cut at 0.8.
    """
    z1, z2 = (float(t) for t in thresholds)
    x = np.asarray(np.ma.filled(np.ma.asarray(value), np.nan))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(x > z1, 0.0,
                        np.where(x > z2, (x - z1) / (z2 - z1), 1.0))


def continuum_test(value, limits):
    if np.isscalar(limits):
        return value >= limits

    low, high = limits
    return (value >= low) & (value <= high)


#========= Continuum-removed feature container ================================

@dataclass
class ContinuumFeature:
    """
    A continuum-removed feature over the native fit range cl1..cr2.

    interior marks the channels between the continuum windows
    (cl2+1..cr1-1): the band minimum/maximum search and the feature
    weight area use only these.
    """
    continuum_type: Literal["linear", "curved"]

    wavelengths: np.ndarray
    continuum_removed: np.ndarray
    continuum: np.ndarray

    left_continuum: np.ndarray
    right_continuum: np.ndarray

    interior: np.ndarray = field(default=None)
    channels: tuple = field(default=None)


#========= specpr / bandmp helpers ============================================

def _wtochbin(wavelengths, w1, w2):
    """
    specpr wtochbin in Tetracorder mode (ier=12): wavelength interval to a
    0-based channel range (ich1, ich2) on the full wavelength array.

    ich1 is the first channel with wavelength >= w1; ich2 the last channel
    from there on with wavelength <= w2. If that is a single channel the
    interval does not actually contain, the channel nearest w1 is used.
    Raises ValueError when w1 is beyond the last channel.
    """
    n = wavelengths.size
    for i in range(n):
        if w1 <= wavelengths[i]:
            ich1 = ich2 = i
            for j in range(i, n):
                if w2 < wavelengths[j]:
                    break
                ich2 = j
            break
    else:
        raise ValueError(f"Wavelength window {w1:g} - {w2:g} is out of range")

    if ich1 == ich2 and not (w1 <= wavelengths[ich1] and w2 >= wavelengths[ich2]):
        if ich1 > 0 and abs(wavelengths[ich1] - w1) > abs(wavelengths[ich1 - 1] - w1):
            return ich1 - 1, ich1 - 1
        if ich1 < n - 1 and abs(wavelengths[ich1] - w1) > abs(wavelengths[ich1 + 1] - w1):
            return ich1 + 1, ich1 + 1
    return ich1, ich2


def _prepare(spectra, wavelengths, valid_bands):
    spectra = np.ma.filled(np.ma.asarray(spectra, dtype=np.float32), np.nan)
    spectra = np.asarray(spectra, dtype=np.float32)
    wavelengths = np.asarray(wavelengths, dtype=np.float32)

    if wavelengths.ndim != 1:
        raise ValueError(f"Wavelengths must be one-dimensional; got shape {wavelengths.shape}.")
    if spectra.shape[-1] != wavelengths.size:
        raise ValueError(
            f"Spectra have {spectra.shape[-1]} bands but "
            f"{wavelengths.size} wavelengths were supplied."
        )

    if valid_bands is None:
        valid_bands = np.ones(wavelengths.shape, dtype=bool)
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)
        if valid_bands.shape != wavelengths.shape:
            raise ValueError("valid_bands must have the same shape as wavelengths")

    return spectra, wavelengths, valid_bands


# ======Refactored for perfomance======================================================
# Slower than the below.This processes every channel of the feature span
# while number of channels per window is almost always significantly less
# than the size of the span provided when called by linear continuum
# neglibible on curved continuum, which pass pre-sliced windows.

# Preserved as summation order differs (NumPy pairwise vs
# channel-by-channel), giving up to ~4 ulp differences in window means. I may revert.

# def _window_mean(values, wav, mask):
#     """
#     REAL*4 means of reflectance and wavelength over mask (last axis),
#     skipping NaN, as the bandmp continuum loops do. Returns means and counts.
#     """
#     ok = mask & np.isfinite(values)
#     n = ok.sum(axis=-1)
#     with np.errstate(invalid="ignore", divide="ignore"):
#         refl = (np.where(ok, values, np.float32(0)).sum(axis=-1, dtype=np.float32)
#                 / n.astype(np.float32))
#         wl = (np.where(ok, wav, np.float32(0)).sum(axis=-1, dtype=np.float32)
#               / n.astype(np.float32))
#     return refl.astype(np.float32), wl.astype(np.float32), n
# =============================================================================

def _window_mean(values, wav, mask):
    """
    REAL*4 means of reflectance and wavelength over mask (last axis),
    skipping NaN, as the bandmp continuum loops do. Returns means and counts.

    Accumulates channel by channel in channel order (the Fortran DO loop),
    touching only the masked channels.
    """
    shape = values.shape[:-1]
    refl_sum = np.zeros(shape, dtype=np.float32)
    wl_sum = np.zeros(shape, dtype=np.float32)
    n = np.zeros(shape, dtype=np.intp)

    for c in np.flatnonzero(mask):
        v = values[..., c]
        good = np.isfinite(v)
        refl_sum += np.where(good, v, np.float32(0))
        wl_sum += np.where(good, np.float32(wav[c]), np.float32(0))
        n += good

    with np.errstate(invalid="ignore", divide="ignore"):
        count = n.astype(np.float32)
        refl = refl_sum / count
        wl = wl_sum / count
    return refl, wl, n

def _remove_continuum(values, continuum, use):
    """rfobsc = rfobs / contin where valid and |contin| > 0.1e-20, else deleted."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(
            use & np.isfinite(values) & (np.abs(continuum) > _TINY),
            values / continuum,
            np.nan,
        ).astype(np.float32)


#========= Continuum removal ==================================================

def linear_feature_continuum(
    spectra,
    wavelengths,
    left_window,
    right_window,
    valid_bands=None,
):
    """
    Linear continuum removal: bdmset for a 1-D reference, bandmp for
    target spectra (..., bands).

    Channels are cl1..cl2 and cr1..cr2 from wtochbin. The continuum line
    goes through the (reflectance, wavelength) means of the two windows,
    each skipping deleted channels. The continuum-removed spectrum covers
    cl1..cr2 (imgflg = 0).

    Where a spectrum has no usable channel in a window, a window mean
    below 0.1e-20, or a zero wavelength span, bandmp deletes its output:
    here its continuum-removed values are NaN.

    Raises ValueError if the windows fall outside the wavelength range or
    are not in sequence (getifeat disables such a feature).
    """
    spectra, wavelengths, valid_bands = _prepare(spectra, wavelengths, valid_bands)

    cl1, cl2 = _wtochbin(wavelengths, float(left_window[0]), float(left_window[1]))
    cr1, cr2 = _wtochbin(wavelengths, float(right_window[0]), float(right_window[1]))
    if cl1 > cl2 or cl2 + 1 >= cr1 or cr1 > cr2:
        raise ValueError("Continuum channels are not in a valid sequence")

    span = slice(cl1, cr2 + 1)
    wav = wavelengths[span]
    k = np.arange(cl1, cr2 + 1)
    left = k <= cl2
    right = k >= cr1
    interior = ~left & ~right

    values = np.where(valid_bands[span], spectra[..., span], np.nan).astype(np.float32)

    left_refl, left_wl, n_left = _window_mean(values, wav, left)
    right_refl, right_wl, n_right = _window_mean(values, wav, right)

    bottom = (right_wl - left_wl).astype(np.float32)
    deleted = ((n_left < 1) | (n_right < 1)
               | ~(left_refl >= _TINY) | ~(right_refl >= _TINY)
               | ~(np.abs(bottom) >= _TINY))

    with np.errstate(invalid="ignore", divide="ignore"):
        slope = ((right_refl - left_refl) / bottom).astype(np.float32)
        intercept = (right_refl - slope * right_wl).astype(np.float32)
        continuum = (slope[..., None] * wav + intercept[..., None]).astype(np.float32)

    continuum_removed = _remove_continuum(values, continuum,
                                          np.ones(wav.shape, dtype=bool))
    continuum_removed[deleted] = np.nan

    return ContinuumFeature(
        continuum_type="linear",
        wavelengths=wav,
        continuum_removed=np.ma.masked_invalid(continuum_removed),
        continuum=continuum,
        left_continuum=left_refl,
        right_continuum=right_refl,
        interior=interior,
        channels=(cl1, cl2, cr1, cr2),
    )


def curved_feature_continuum(spectra, wavelengths, continuum_windows, valid_bands=None):
    """
    Curved continuum removal for target spectra: bandmpcv.

    continuum_windows:
        [
            (outer_left_start, outer_left_stop),
            (inner_left_start, inner_left_stop),
            (inner_right_start, inner_right_stop),
            (outer_right_start, outer_right_stop),
        ]

    A natural cubic spline (icsicu with zero end conditions) goes through
    the four window means at the window mid-wavelengths (getifeat cvwav).
    The continuum-removed spectrum covers the inner windows and the band,
    inner_left_start..inner_right_stop (imgflg = 0).

    valid_bands is honoured in the two inner windows only. bandmpcv
    skips only deleted OBSERVED points in the outer windows (it tests
    rflibc there, but rlbc is 0.0 outside the inner windows,
    getifeat.r:1674-1710), so [DELETPTS] channels are averaged into
    the outer anchors.

    NOTE: natively the REFERENCE of a curved feature is continuum removed
    by bdmset with a straight line through the two inner windows. Use
    linear_feature_continuum(reference, wavelengths, windows[1], windows[2])
    for the reference; it gives the same channel range.
    """
    spectra, wavelengths, valid_bands = _prepare(spectra, wavelengths, valid_bands)

    if len(continuum_windows) != 4:
        raise ValueError("Curved continuum requires four wavelength windows")

    ch = [_wtochbin(wavelengths, float(a), float(b)) for a, b in continuum_windows]
    (n1, n2), (n3, n4), (n5, n6), (n7, n8) = ch
    if (n1 > n2 or n2 >= n3 or n3 > n4 or n4 + 1 >= n5
            or n5 > n6 or n6 >= n7 or n7 > n8):
        raise ValueError("Curved continuum channels are not in a valid sequence")

    cwav = np.array([(float(a) + float(b)) / 2.0 for a, b in continuum_windows],
                    dtype=np.float32)
    if not np.all(np.diff(cwav) > 0):
        raise ValueError("Curved continuum window mid-wavelengths are not consecutive")

    means = []
    counts = []
    wl_means = []
    for i, (a, b) in enumerate(ch):
        idx = slice(a, b + 1)
        vals = spectra[..., idx]
        if i in (1, 2):
            vals = np.where(valid_bands[idx], vals, np.nan).astype(np.float32)
        m, w, n = _window_mean(vals, wavelengths[idx], np.ones(b - a + 1, dtype=bool))
        means.append(m)
        counts.append(n)
        wl_means.append(w)

    anchors = np.stack(means, axis=-1).astype(np.float32)          # (..., 4)
    deleted = np.zeros(anchors.shape[:-1], dtype=bool)
    for m, n in zip(means, counts):
        deleted |= (n < 1) | ~(m >= _TINY)
    # the linear-continuum wavelength span check is still done in bandmpcv
    deleted |= ~(np.abs(wl_means[2] - wl_means[1]) >= _TINY)

    span = slice(n3, n6 + 1)
    wav = wavelengths[span]
    k = np.arange(n3, n6 + 1)
    interior = (k > n4) & (k < n5)

    safe = np.where(deleted[..., None], np.float32(1.0), anchors)
    spline = CubicSpline(cwav.astype(float), safe.astype(float), axis=-1,
                         bc_type="natural")
    continuum = np.asarray(spline(wav.astype(float)), dtype=np.float32)

    values = np.where(valid_bands[span], spectra[..., span], np.nan).astype(np.float32)
    continuum_removed = _remove_continuum(values, continuum,
                                          np.ones(wav.shape, dtype=bool))
    continuum_removed[deleted] = np.nan

    return ContinuumFeature(
        continuum_type="curved",
        wavelengths=wav,
        continuum_removed=np.ma.masked_invalid(continuum_removed),
        continuum=continuum,
        left_continuum=means[1],
        right_continuum=means[2],
        interior=interior,
        channels=(n3, n4, n5, n6),
    )


#========= Characterise primitive =============================================

def characterise_feature(
    reference,
    target,
    polarity="absorption",
):
    """
    Least-squares fit of a continuum-removed reference to continuum-removed
    target spectra: bdmset band min/max + bandmp/bandmpcv sums, plus the
    per-feature checks at the top of tp1mat.

    Parameters
    ----------
    reference : ContinuumFeature
        1-D reference (linear_feature_continuum).
    target : ContinuumFeature
        Target spectra (..., bands), same channel range as the reference.
    polarity : {"absorption", "emission", None}
        Accepted for compatibility. As in bdmset, the feature type is
        always taken from the reference: emission when (rmax - 1) >
        (1 - rmin) over the band channels, otherwise absorption.

    Returns
    -------
    dict
        fit and depth (masked where bandmp deletes its output or tp1mat
        zeroes the feature), continuum at the band extremum, left/right
        continuum and continuum-shape ratios.

    Masked (feature zeroed) where:
        bandmp:  no channels summed, |b| < 0.1e-20, |xk + 1| < 0.1e-20,
                 or a deleted continuum;
        tp1mat:  fit outside (0, 1.1), |depth| < 0.1e-6, or continuum,
                 left or right continuum outside the default ct/lct/rct
                 limits 0.1e-6 .. 0.1e20.
    """
    if polarity not in {None, "absorption", "emission"}:
        raise ValueError("Polarity must be None, 'absorption' or 'emission'")

    reference_wl = np.asarray(reference.wavelengths, dtype=np.float32)
    target_wl = np.asarray(target.wavelengths, dtype=np.float32)
    if reference_wl.shape != target_wl.shape or not np.allclose(reference_wl, target_wl):
        raise ValueError("Reference and target feature wavelengths do not match.")

    rflibc = np.asarray(np.ma.filled(np.ma.asarray(reference.continuum_removed,
                                                   dtype=np.float32), np.nan))
    rfobsc = np.asarray(np.ma.filled(np.ma.asarray(target.continuum_removed,
                                                   dtype=np.float32), np.nan))
    if rflibc.ndim != 1:
        raise ValueError("Reference feature must be one-dimensional.")
    if rfobsc.shape[-1] != rflibc.size:
        raise ValueError("Reference and target features contain different numbers of bands.")

    interior = reference.interior
    if interior is None:
        interior = np.ones(rflibc.shape, dtype=bool)

    # ------------------------------------------------------------------
    # bdmset: band minimum / maximum over the band channels only,
    # first occurrence of a strictly smaller / larger value.
    # ------------------------------------------------------------------
    inner = np.flatnonzero(interior & np.isfinite(rflibc))
    if inner.size == 0:
        raise ValueError("All reference band channels are deleted.")

    minch = inner[np.argmin(rflibc[inner])]
    maxch = inner[np.argmax(rflibc[inner])]
    rmin = rflibc[minch]
    rmax = rflibc[maxch]
    if rmin < 0.0:
        raise ValueError("Reference band depth is less than zero.")

    if (rmax - 1.0) > (1.0 - rmin):
        polarity, extremum_index = "emission", int(maxch)
    else:
        polarity, extremum_index = "absorption", int(minch)

    reference_depth = np.float32(1.0 - rflibc[extremum_index])

    # ------------------------------------------------------------------
    # bandmp sums over cl1..cr2, REAL*8
    # ------------------------------------------------------------------

    # There was a significant refactor here to avoid broadcasting to p x n
    # intermediate arrays and take advantage of matrix-vector products
    # original not preserved as cast to float32 should eliminate any
    # difference in the float64 values from accumulation order.
    ref_ok = np.isfinite(rflibc)
    r = np.where(ref_ok, rflibc, 0).astype(np.float64)       # 1-D, REAL*8
    ok = np.isfinite(rfobsc) & ref_ok                        # P * n bool
    n = ok.sum(axis=-1)
    okf = ok.astype(np.float64)                              
    dro = np.where(ok, rfobsc, 0).astype(np.float64)         

    suml = okf @ r
    sumll = okf @ (r * r)
    sumol = dro @ r
    sumo = dro.sum(axis=-1)
    sumoo = np.einsum("...i,...i->...", dro, dro)
    dxn = np.where(n > 0, n, 1).astype(np.float64)



    top = (sumol - sumo * suml / dxn).astype(np.float32)
    bottom = (sumll - suml * suml / dxn).astype(np.float32)
    botm2 = (sumoo - sumo * sumo / dxn).astype(np.float32)

    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        slope = np.where(np.abs(bottom) < _TINY, 0, top / bottom).astype(np.float32)
        xk = ((1.0 - slope) / slope).astype(np.float32)
        xk1 = (xk + 1.0).astype(np.float32)
        rftemp = ((rflibc[extremum_index] + xk) / xk1).astype(np.float32)
        depth = (1.0 - rftemp).astype(np.float32)
        bprime = np.where(np.abs(botm2) < _TINY, 0, top / botm2).astype(np.float32)
        fit = np.sqrt(np.abs(slope * bprime)).astype(np.float32)

    # ------------------------------------------------------------------
    # Target continuum at the band extremum (conref) and the window means
    # ------------------------------------------------------------------
    continuum = np.asarray(target.continuum, dtype=np.float32)[..., extremum_index]
    left_continuum = np.asarray(target.left_continuum, dtype=np.float32)
    right_continuum = np.asarray(target.right_continuum, dtype=np.float32)

    deleted = (
        (n == 0)
        | ~(np.abs(slope) >= _TINY)
        | ~(np.abs(xk1) >= _TINY)
        # tp1mat, before any rule test
        | ~((fit > 0.0) & (fit < 1.1))
        | ~(np.abs(depth) >= 0.1e-6)
    )
    for c in (continuum, left_continuum, right_continuum):
        deleted |= ~((c >= _CONT_MIN) & (c <= _CONT_MAX))

    fit = np.ma.masked_where(deleted, fit)
    depth = np.ma.masked_where(deleted, depth)
    fit_depth = np.ma.asarray(fit * depth, dtype=np.float32)

    band_bottom = np.ma.asarray(continuum * (np.float32(1.0) - depth), dtype=np.float32)
    reflectance_depth = np.ma.asarray(continuum * depth, dtype=np.float32)

    # ------------------------------------------------------------------
    # Continuum-shape ratios with the tp1mat clamps
    #   lct/rct, rct/lct:   denominator >= 0.1e-6, numerator <= 0.1e20
    #   shoulder ratios:    |denominator| < 0.1e-6 -> +0.1e-6,
    #                       numerator <= 0.1e20
    # ------------------------------------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        left_right_ratio = (np.minimum(left_continuum, _CONT_MAX)
                            / np.maximum(right_continuum, _CONT_MIN))
        right_left_ratio = (np.minimum(right_continuum, _CONT_MAX)
                            / np.maximum(left_continuum, _CONT_MIN))

        def shoulder(numerator_side, denominator_side):
            den = denominator_side - band_bottom
            den = np.ma.where(np.abs(den) < _CONT_MIN, _CONT_MIN, den)
            return np.minimum(numerator_side - band_bottom, _CONT_MAX) / den

        left_shoulder_ratio = shoulder(left_continuum, right_continuum)
        right_shoulder_ratio = shoulder(right_continuum, left_continuum)

    return {
        "fit": fit,
        "depth": depth,
        "fit_depth": fit_depth,
        "slope": slope,

        "polarity": polarity,
        "extremum_index": extremum_index,
        "reference_depth": reference_depth,

        "continuum": continuum,
        "left_continuum": left_continuum,
        "right_continuum": right_continuum,
        "band_bottom": band_bottom,
        "reflectance_depth": reflectance_depth,

        "left_right_ratio": left_right_ratio,
        "right_left_ratio": right_left_ratio,
        "left_shoulder_ratio": left_shoulder_ratio,
        "right_shoulder_ratio": right_shoulder_ratio,

        # getifeat weight area: band channels only
        "reference_area": float(np.nansum(np.abs(1.0 - rflibc[interior]))),
        "value_count": int(np.isfinite(rflibc).sum()),
    }
