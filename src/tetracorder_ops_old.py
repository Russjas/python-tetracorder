"""tetracorder operations reimplemented in numpy scipy.
to be used by the rule interpreter algorithm once implemented"""

from dataclasses import dataclass
from typing import Literal
from copy import deepcopy
from pathlib import Path
import re
import sqlite3

import numpy as np
from scipy.interpolate import CubicSpline

from .config import VARIABLE_PRESETS

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


def load_references(db_file: str|Path) ->dict:
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
            Spectra.YData
        FROM Samples
        JOIN Spectra USING (SampleID)
        WHERE Samples.ConvolvedRecord IS NOT NULL
        """
    )

    for row in rows:

        key = (
            row["Library"],
            row["ConvolvedRecord"],
        )

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
            "wavelengths": np.frombuffer(
                row["XData"],
                dtype=np.float32,
            ).copy(),
            "reflectance": np.frombuffer(
                row["YData"],
                dtype=np.float32,
            ).copy(),
        }

    connection.close()

    return references


#========= Helpers for the rule evaluation logic==============================

def fuzzy_greater(value, thresholds):
    reject, full = thresholds
    value = np.ma.asarray(value)
    result = np.zeros(value.shape, dtype=float)

    valid = np.isfinite(np.ma.filled(value, np.nan))

    result[valid & (value >= full)] = 1.0

    transition = (valid & (value > reject) & (value < full))

    result[transition] = ((value[transition] - reject) / (full - reject))

    return result

def fuzzy_less(value, thresholds):
    reject, full = thresholds

    value = np.ma.asarray(value)
    result = np.zeros(value.shape, dtype=float)

    valid = np.isfinite(np.ma.filled(value, np.nan))

    result[valid & (value <= full)] = 1.0

    transition = (valid & (value < reject) & (value > full))

    result[transition] = ((reject - value[transition]) / (reject - full))

    return result

def continuum_test(value, limits):
    if np.isscalar(limits):
        return value >= limits

    low, high = limits
    return (value >= low) & (value <= high)

#==========================================================================

class GaussianConvolver:
    """
    Convolve high-resolution library spectra to scanner channels using
    Gaussian spectral response functions defined by scanner FWHM.

    The convolution weights are calculated once and can then be reused
    for multiple library spectra on the same wavelength grid.

    Usage: 
    convolver = GaussianConvolver(
                lib_wl=library_wavelengths,
                scanner_wl=scanner_wavelengths,
                scanner_fwhm=scanner_fwhm,
                )
    convolved_reflectance = convolver.convolve(library_reflectance,)
    """

    def __init__(self, lib_wl, scanner_wl, scanner_fwhm):
        self.lib_wl = np.asarray(lib_wl, dtype=float)
        self.scanner_wl = np.asarray(scanner_wl, dtype=float)
        self.scanner_fwhm = np.asarray(scanner_fwhm, dtype=float)

        if self.scanner_fwhm.ndim == 0 or self.scanner_fwhm.size == 1:
            self.scanner_fwhm = np.full(
                self.scanner_wl.shape,
                float(self.scanner_fwhm.ravel()[0]),
            )

        if self.scanner_wl.shape != self.scanner_fwhm.shape:
            raise ValueError(
                "scanner_wl and scanner_fwhm must have equal lengths"
            )

        if self.lib_wl.ndim != 1:
            raise ValueError("lib_wl must be one-dimensional")

        if np.any(np.diff(self.lib_wl) <= 0):
            raise ValueError("lib_wl must be strictly increasing")

        if np.any(self.scanner_fwhm <= 0):
            raise ValueError("All scanner FWHM values must be positive")

        self.weights = self._build_weights()

    def _build_weights(self):
        # Gaussian response:
        # exp(-4 ln(2) * ((x - centre) / FWHM)**2)
        offset = (
            self.lib_wl[None, :]
            - self.scanner_wl[:, None]
        )

        weights = np.exp(
            -4.0 * np.log(2.0)
            * (offset / self.scanner_fwhm[:, None]) ** 2
        )

        # Approximate integration widths on the library wavelength grid.
        # This matters if the library grid is not perfectly uniform.
        spacing = np.gradient(self.lib_wl)
        weights *= spacing[None, :]

        # Do not calculate scanner bands whose centres fall outside the
        # library wavelength coverage.
        outside = (
            (self.scanner_wl <= self.lib_wl.min())
            | (self.scanner_wl >= self.lib_wl.max())
        )
        weights[outside] = 0.0

        # Negligible weights add computation without affecting the result.
        weights[weights < 1e-8] = 0.0

        return weights

    def convolve(self, lib_refl):
        """
        Parameters
        ----------
        lib_refl : ndarray
            Either:
                (library_bands,)
            or:
                (number_of_spectra, library_bands)

        Returns
        -------
        ndarray
            Spectrum/spectra convolved to scanner_wl.
        """
        spectra = np.asarray(lib_refl, dtype=float)
        was_1d = spectra.ndim == 1

        if was_1d:
            spectra = spectra[None, :]

        if spectra.ndim != 2:
            raise ValueError("lib_refl must be one- or two-dimensional")

        if spectra.shape[1] != self.lib_wl.size:
            raise ValueError(
                f"Expected {self.lib_wl.size} library bands, "
                f"received {spectra.shape[1]}"
            )

        valid = np.isfinite(spectra)

        # Renormalise independently where individual spectra contain NaNs.
        numerator = (
            np.where(valid, spectra, 0.0)
            @ self.weights.T
        )

        denominator = (
            valid.astype(float)
            @ self.weights.T
        )

        result = np.full(numerator.shape, np.nan, dtype=float)

        np.divide(
            numerator,
            denominator,
            out=result,
            where=denominator > 0,
        )

        if was_1d:
            return result[0]

        return result


@dataclass
class ContinuumFeature:
    """A class for holding a continuum removed feature for the python tetracorder algorithm"""
    continuum_type: Literal["linear", "curved"]

    wavelengths: np.ndarray
    continuum_removed: np.ndarray
    continuum: np.ndarray

    left_continuum: np.ndarray
    right_continuum: np.ndarray
# Validation for wavelenth dataset

def _validate_wavelengths(wavelengths):
    """
    Validate and return a one-dimensional, increasing wavelength array.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)

    if wavelengths.ndim != 1:
        raise ValueError(
            f"Wavelengths must be one-dimensional; got shape "
            f"{wavelengths.shape}."
        )

    if wavelengths.size == 0:
        raise ValueError("Wavelength array is empty.")

    if not np.all(np.isfinite(wavelengths)):
        raise ValueError("Wavelength array contains non-finite values.")

    if not np.all(np.diff(wavelengths) > 0):
        raise ValueError(
            "Wavelengths must be strictly increasing with no duplicates."
        )

    return wavelengths

def _validate_wavelengths(wavelengths, valid_bands=None):
    """
    Validate wavelength data.

    The full wavelength array may contain non-monotonic channels if those
    channels are excluded by valid_bands. The usable wavelength sequence
    must be strictly increasing.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)

    if wavelengths.ndim != 1:
        raise ValueError(f"Wavelengths must be one-dimensional; got shape{wavelengths.shape}.")

    if wavelengths.size == 0:
        raise ValueError("Wavelength array is empty.")

    if not np.all(np.isfinite(wavelengths)):
        raise ValueError("Wavelength array contains non-finite values.")

    if valid_bands is None:
        selected = wavelengths
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)

        if valid_bands.shape != wavelengths.shape:
            raise ValueError("valid_bands must have the same shape as wavelengths")

        selected = wavelengths[valid_bands]

    if not np.all(np.diff(selected) > 0):
        raise ValueError("Valid wavelengths must be strictly increasing with no duplicates.")

    return wavelengths


def _validate_window(wavelengths, window, name="window", valid_bands = None):
    """
    Check that a window contains two valid bounds and at least one band.

    Returns
    -------
    valid : bool
    message : str
    """
    try:
        lower, upper = window
    except (TypeError, ValueError):
        return (False, f"{name} must contain exactly two wavelength bounds; got {window!r}.",)

    try:
        lower = float(lower)
        upper = float(upper)
    except (TypeError, ValueError):
        return (False, f"{name} bounds must be numeric; got {window!r}.",)

    if not np.isfinite(lower) or not np.isfinite(upper):
        return (False, f"{name} contains a non-finite wavelength bound: {window!r}.",)

    if lower >= upper:
        return (False, f"{name} lower bound must be less than its upper bound; got {lower:g}–{upper:g}.",)

    mask = (wavelengths >= lower) & (wavelengths <= upper)

    if not np.any(mask):
        return (False, f"{name} {lower:g} - {upper:g} contains no bands.",)
    if valid_bands is not None:
        mask &= valid_bands

        if not np.any(mask):
            return (False, f"{name} {lower:g} - {upper:g} contains no valid bands.",)
    return True, ""

def _window_mean(spectra, wavelengths, window, name="window", valid_bands = None):
    lower, upper = window

    selected = (
        (wavelengths >= lower)
        & (wavelengths <= upper)
    )
    if valid_bands is not None:
        selected &= valid_bands
    values = spectra[..., selected]

    finite = np.isfinite(values)
    count = finite.sum(axis=-1)

    total = np.where(finite, values, 0.0).sum(axis=-1)

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count

    mean = np.where(count > 0, mean, np.nan)

    return mean


# Continuum calculation ported from specspr

def linear_feature_continuum(
    spectra,
    wavelengths,
    left_window,
    right_window,
    valid_bands=None,
):
    """
    Reproduce the linear continuum used by Specpr/Tetracorder.

    spectra:
        (..., bands)

    wavelengths:
        (bands,)

    Returns
    -------
    ContinuumFeature
        Continuum-removed feature data and associated continuum geometry.

    Notes
    -----
    Tetracorder bandmp/bdmset operate on REAL*4 spectra and continuum
    arithmetic.

    Deleted/invalid points are omitted from both the reflectance and
    wavelength continuum-window averages, matching the Ratfor logic.
    """

    spectra = np.ma.asarray(
        spectra,
        dtype=np.float32,
    )

    wavelengths = np.asarray(
        wavelengths,
        dtype=np.float32,
    )

    if valid_bands is None:
        valid_bands = np.ones(
            wavelengths.shape,
            dtype=bool,
        )
    else:
        valid_bands = np.asarray(
            valid_bands,
            dtype=bool,
        )

    if valid_bands.shape != wavelengths.shape:
        raise ValueError(
            "valid_bands must have the same shape as wavelengths"
        )

    wavelengths = _validate_wavelengths(
        wavelengths,
        valid_bands=valid_bands,
    )

    # _validate_wavelengths currently converts to ordinary Python
    # float/float64 internally. Restore the REAL*4-equivalent dtype.
    wavelengths = np.asarray(
        wavelengths,
        dtype=np.float32,
    )

    if spectra.shape[-1] != wavelengths.size:
        raise ValueError(
            f"Spectra have {spectra.shape[-1]} bands but "
            f"{wavelengths.size} wavelengths were supplied."
        )

    check, msg = _validate_window(
        wavelengths,
        left_window,
        name="Left continuum window",
        valid_bands=valid_bands,
    )
    if not check:
        raise ValueError(msg)

    check, msg = _validate_window(
        wavelengths,
        right_window,
        name="Right continuum window",
        valid_bands=valid_bands,
    )
    if not check:
        raise ValueError(msg)

    left = (
        (wavelengths >= left_window[0])
        & (wavelengths <= left_window[1])
        & valid_bands
    )

    right = (
        (wavelengths >= right_window[0])
        & (wavelengths <= right_window[1])
        & valid_bands
    )

    # bandmp imaging mode:
    #
    #   il = cl2 + 1
    #   ir = cr1 - 1
    #
    # so continuum-window channels are excluded from the fitted feature.
    feature = (
        (wavelengths > left_window[1])
        & (wavelengths < right_window[0])
        & valid_bands
    )

    if not feature.any():
        raise ValueError(
            "No bands in feature interval"
        )

    # ------------------------------------------------------------------
    # Ratfor-style continuum-window means.
    #
    # bandmp/bdmset:
    #
    #   if (wav(i) == delpt || rfobs(i) == delpt) next
    #
    #   avlc  = avlc  + rfobs(i)
    #   avwlc = avwlc + wav(i)
    #
    # so a deleted reflectance channel is omitted from BOTH the
    # reflectance average and the corresponding wavelength average.
    # ------------------------------------------------------------------

    def continuum_mean(mask):

        values = np.ma.asarray(
            spectra[..., mask],
            dtype=np.float32,
        )

        data = np.asarray(
            np.ma.filled(
                values,
                np.nan,
            ),
            dtype=np.float32,
        )

        window_wl = np.asarray(
            wavelengths[mask],
            dtype=np.float32,
        )

        finite = np.isfinite(data)

        count = finite.sum(
            axis=-1,
        )

        refl_total = np.sum(
            np.where(
                finite,
                data,
                np.float32(0.0),
            ),
            axis=-1,
            dtype=np.float32,
        )

        # window_wl is 1-D and broadcasts over all leading spectrum
        # dimensions, so each target spectrum gets the wavelength mean
        # corresponding to its own surviving continuum channels.
        wl_total = np.sum(
            np.where(
                finite,
                window_wl,
                np.float32(0.0),
            ),
            axis=-1,
            dtype=np.float32,
        )

        refl_mean = np.full(
            np.shape(count),
            np.nan,
            dtype=np.float32,
        )

        wl_mean = np.full(
            np.shape(count),
            np.nan,
            dtype=np.float32,
        )

        count32 = count.astype(
            np.float32,
        )

        np.divide(
            refl_total,
            count32,
            out=refl_mean,
            where=count > 0,
        )

        np.divide(
            wl_total,
            count32,
            out=wl_mean,
            where=count > 0,
        )

        return refl_mean, wl_mean

    left_refl, left_wl = continuum_mean(
        left
    )

    right_refl, right_wl = continuum_mean(
        right
    )

    # ------------------------------------------------------------------
    # Continuum line.
    #
    # Ratfor:
    #
    #   bottom = avwrc - avwlc
    #   a = (avrc - avlc) / bottom
    #   b = avrc - a * avwrc
    #
    # All of these are REAL*4.
    # ------------------------------------------------------------------

    bottom = np.asarray(
        right_wl - left_wl,
        dtype=np.float32,
    )

    # A spectrum with no valid continuum points has NaN bottom and will
    # naturally propagate as invalid. This catches the native
    # effectively-zero wavelength span condition.
    bad_bottom = (
        np.isfinite(bottom)
        & (
            np.abs(bottom)
            < np.float32(1.0e-21)
        )
    )

    if np.any(bad_bottom):
        raise ValueError(
            "Wavelength range of continuum is too small"
        )

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        slope = np.asarray(
            (
                right_refl
                - left_refl
            )
            / bottom,
            dtype=np.float32,
        )

    intercept = np.asarray(
        right_refl
        - slope * right_wl,
        dtype=np.float32,
    )

    feature_wl = wavelengths[
        feature
    ]

    continuum = np.asarray(
        slope[..., None]
        * feature_wl
        + intercept[..., None],
        dtype=np.float32,
    )

    # ------------------------------------------------------------------
    # Continuum removal.
    #
    # Ratfor:
    #
    #   contin = a * wav(i) + b
    #
    #   if (abs(contin) > 0.1e-20)
    #       rfobsc(i) = rfobs(i) / contin
    #   else
    #       rfobsc(i) = delpt
    # ------------------------------------------------------------------

    feature_values = np.ma.asarray(
        spectra[..., feature],
        dtype=np.float32,
    )

    feature_data = np.asarray(
        np.ma.filled(
            feature_values,
            np.nan,
        ),
        dtype=np.float32,
    )

    continuum_removed = np.full(
        continuum.shape,
        np.nan,
        dtype=np.float32,
    )

    valid = (
        np.isfinite(feature_data)
        & np.isfinite(continuum)
        & (
            np.abs(continuum)
            >= np.float32(1.0e-21)
        )
    )

    np.divide(
        feature_data,
        continuum,
        out=continuum_removed,
        where=valid,
    )

    continuum_removed = np.ma.masked_invalid(
        continuum_removed
    )

    return ContinuumFeature(
        continuum_type="linear",
        wavelengths=feature_wl,
        continuum_removed=continuum_removed,
        continuum=continuum,
        left_continuum=left_refl,
        right_continuum=right_refl,
    )


def curved_feature_continuum(spectra, wavelengths, continuum_windows, valid_bands = None):
    """
    Reproduce the curved continuum used by Specpr/Tetracorder.

    continuum_windows:
        [
            (outer_left_start, outer_left_stop),
            (inner_left_start, inner_left_stop),
            (inner_right_start, inner_right_stop),
            (outer_right_start, outer_right_stop),
        ]
    """
    spectra = np.asanyarray(spectra)
    

    if valid_bands is None:
        valid_bands = np.ones(wavelengths.shape, dtype=bool)
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)

    if valid_bands.shape != wavelengths.shape:
        raise ValueError("valid_bands must have the same shape as wavelengths")

    wavelengths = _validate_wavelengths(wavelengths, valid_bands=valid_bands)

    if spectra.shape[-1] != wavelengths.size:
        raise ValueError(
            f"Spectra have {spectra.shape[-1]} bands but "
            f"{wavelengths.size} wavelengths were supplied."
        )

    if len(continuum_windows) != 4:
        raise ValueError(
            "Curved continuum requires four wavelength windows"
        )

    for i, window in enumerate(continuum_windows, start=1):
        check, msg = _validate_window(
            wavelengths, window,
            name=f"Continuum window {i}", valid_bands = valid_bands)
        if not check:
            raise ValueError(msg)


    anchor_wavelengths = np.array([(start + stop) / 2.0
            for start, stop in continuum_windows
        ],
        dtype=float,
    )

    anchor_reflectances = []

    for i, window in enumerate(continuum_windows, 1):
        anchor = _window_mean(
            spectra, wavelengths,
            window, name=f"Continuum window {i}",
            valid_bands = valid_bands)
    
        anchor_reflectances.append(anchor)
    
    anchor_reflectances = np.stack(anchor_reflectances,axis=-1)
    
    
    if not np.all(np.isfinite(anchor_reflectances)):
        raise ValueError(
            "One or more continuum windows contain no finite reflectance values")
        
    feature = ((wavelengths > continuum_windows[1][1]) & (wavelengths < continuum_windows[2][0])
               & valid_bands)

    if not feature.any():
        raise ValueError("No bands in feature interval")

    feature_wl = wavelengths[feature]

    spline = CubicSpline(
        anchor_wavelengths,
        anchor_reflectances,
        axis=-1,
        bc_type="natural",
    )

    continuum = spline(feature_wl)
    
    with np.errstate(divide="ignore", invalid="ignore"):
        continuum_removed = spectra[..., feature] / continuum

    continuum_removed = np.ma.masked_invalid(continuum_removed)

    return ContinuumFeature(
    continuum_type="curved",
    wavelengths=feature_wl,
    continuum_removed=continuum_removed,
    continuum=continuum,
    left_continuum=anchor_reflectances[..., 1],
    right_continuum=anchor_reflectances[..., 2],
)

# ===============Characterise primitive ========================
def characterise_feature(
    reference,
    target,
    polarity="absorption",
):
    """
    Characterise the relationship between a reference feature and one
    or more target features using the Specpr/Tetracorder fitting method.

    Parameters
    ----------
    reference : ContinuumFeature
        A single reference feature. Its continuum_removed array must
        have shape (bands,).

    target : ContinuumFeature
        Target feature or features. Its continuum_removed array may
        have shape (bands,) or (..., bands).

    polarity : {"absorption", "emission", None}
        Expected feature polarity. If None, infer polarity from the
        reference as bdmset does.

    Returns
    -------
    dict
        Feature fit, depth and continuum-shape measurements.

    Notes
    -----
    This reproduces the numerical structure of bandmp:

        rflibc, rfobsc       REAL*4
        suml, sumll,
        sumol, sumo, sumoo   REAL*8
        top, bottom, botm2   REAL*4
        b, bprime, fit,
        depth                REAL*4

    Target-invalid channels are also excluded from the corresponding
    reference sums for that target spectrum, matching the bandmp loop.
    """

    if polarity not in {
        None,
        "absorption",
        "emission",
    }:
        raise ValueError(
            "Polarity must be None, 'absorption' or 'emission'"
        )

    if reference.continuum_type != target.continuum_type:
        raise ValueError(
            "Reference and target use different continuum types"
        )

    reference_wl = np.asarray(
        reference.wavelengths,
        dtype=np.float32,
    )

    target_wl = np.asarray(
        target.wavelengths,
        dtype=np.float32,
    )

    if (
        reference_wl.shape != target_wl.shape
        or not np.allclose(
            reference_wl,
            target_wl,
        )
    ):
        raise ValueError(
            "Reference and target feature wavelengths do not match."
        )

    reference_cr = np.ma.asarray(
        reference.continuum_removed,
        dtype=np.float32,
    )

    target_cr = np.ma.asarray(
        target.continuum_removed,
        dtype=np.float32,
    )

    if reference_cr.ndim != 1:
        raise ValueError(
            "Reference feature must be one-dimensional."
        )

    if target_cr.shape[-1] != reference_cr.size:
        raise ValueError(
            "Reference and target features contain "
            "different numbers of bands."
        )

    if reference_cr.count() < 2:
        raise ValueError(
            "Reference feature has fewer than two valid bands."
        )

    # ------------------------------------------------------------------
    # Determine feature polarity.
    #
    # bdmset:
    #
    #   depth = 1 - rmin
    #   emiss = rmax - 1
    #
    # and emission is selected only when emiss > depth.
    # ------------------------------------------------------------------

    if polarity is None:
        max_absorption = np.float32(
            1.0 - float(np.ma.min(reference_cr))
        )

        max_emission = np.float32(
            float(np.ma.max(reference_cr)) - 1.0
        )

        polarity = (
            "emission"
            if max_emission > max_absorption
            else "absorption"
        )

    if polarity == "absorption":
        extremum_index = int(
            np.ma.argmin(reference_cr)
        )
    else:
        extremum_index = int(
            np.ma.argmax(reference_cr)
        )

    # This is the bdmset minch/maxch reference depth.
    reference_depth = np.float32(
        1.0
        - np.float32(
            reference_cr[extremum_index]
        )
    )

    # ------------------------------------------------------------------
    # Convert deleted/masked points to NaN for vectorised validity
    # handling.
    # ------------------------------------------------------------------

    reference_data = np.asarray(
        np.ma.filled(
            reference_cr,
            np.nan,
        ),
        dtype=np.float32,
    )

    target_data = np.asarray(
        np.ma.filled(
            target_cr,
            np.nan,
        ),
        dtype=np.float32,
    )

    # Equivalent to:
    #
    #   if (rfobsc(i) == delpt) next
    #
    # rflibc must also be valid. Because reference_data is 1-D,
    # NumPy broadcasts it over every target spectrum.
    valid = (
        np.isfinite(target_data)
        & np.isfinite(reference_data)
    )

    n = valid.sum(
        axis=-1,
    )

    ref = np.where(
        valid,
        reference_data,
        np.float32(0.0),
    )

    obs = np.where(
        valid,
        target_data,
        np.float32(0.0),
    )

    # ------------------------------------------------------------------
    # bandmp sums.
    #
    # Ratfor explicitly converts rflibc/rfobsc to DOUBLE before
    # accumulation.
    # ------------------------------------------------------------------

    ref64 = ref.astype(
        np.float64,
        copy=False,
    )

    obs64 = obs.astype(
        np.float64,
        copy=False,
    )

    suml = np.sum(
        ref64,
        axis=-1,
        dtype=np.float64,
    )

    sumll = np.sum(
        ref64 * ref64,
        axis=-1,
        dtype=np.float64,
    )

    sumol = np.sum(
        obs64 * ref64,
        axis=-1,
        dtype=np.float64,
    )

    sumo = np.sum(
        obs64,
        axis=-1,
        dtype=np.float64,
    )

    sumoo = np.sum(
        obs64 * obs64,
        axis=-1,
        dtype=np.float64,
    )

    dxn = n.astype(
        np.float64,
    )

    # Prevent the vectorised calculation itself from dividing by zero.
    # n == 0 is separately marked invalid below.
    safe_n = np.where(
        n > 0,
        dxn,
        1.0,
    )

    # ------------------------------------------------------------------
    # Ratfor:
    #
    #   top    = sngl(sumol - sumo*suml/dxn)
    #   bottom = sngl(sumll - suml*suml/dxn)
    #
    # Both deliberately return to REAL*4 here.
    # ------------------------------------------------------------------

    top = np.asarray(
        sumol
        - (
            sumo
            * suml
            / safe_n
        ),
        dtype=np.float32,
    )

    bottom = np.asarray(
        sumll
        - (
            suml
            * suml
            / safe_n
        ),
        dtype=np.float32,
    )

    # ------------------------------------------------------------------
    # b = forward least-squares slope.
    #
    # Ratfor:
    #
    #   if abs(bottom) < 0.1e-20
    #       b = 0
    #   else
    #       b = top / bottom
    # ------------------------------------------------------------------

    slope = np.zeros(
        np.shape(top),
        dtype=np.float32,
    )

    np.divide(
        top,
        bottom,
        out=slope,
        where=(
            np.abs(bottom)
            >= np.float32(1.0e-21)
        ),
    )

    invalid_slope = (
        (n == 0)
        | (
            np.abs(slope)
            < np.float32(1.0e-21)
        )
    )

    # ------------------------------------------------------------------
    # Reverse fit slope.
    #
    # Ratfor:
    #
    #   botm2 = sngl(sumoo - sumo*sumo/dxn)
    #
    #   if abs(botm2) < 0.1e-20
    #       bprime = 0
    #   else
    #       bprime = top / botm2
    # ------------------------------------------------------------------

    botm2 = np.asarray(
        sumoo
        - (
            sumo
            * sumo
            / safe_n
        ),
        dtype=np.float32,
    )

    reverse_slope = np.zeros(
        np.shape(top),
        dtype=np.float32,
    )

    np.divide(
        top,
        botm2,
        out=reverse_slope,
        where=(
            np.abs(botm2)
            >= np.float32(1.0e-21)
        ),
    )

    # ------------------------------------------------------------------
    # Native goodness of fit:
    #
    #   rfit2 = abs(bb * bprime)
    #   rfit  = sqrt(rfit2)
    #
    # No explicit clipping to 0..1 in bandmp.
    # ------------------------------------------------------------------

    fit = np.asarray(
        np.sqrt(
            np.abs(
                slope
                * reverse_slope
            )
        ),
        dtype=np.float32,
    )

    # ------------------------------------------------------------------
    # Native band depth.
    #
    # bandmp constructs:
    #
    #   xk = (1-b)/b
    #   rftemp = (rflibc+xk)/(xk+1)
    #   bd = 1-rftemp(minch)
    #
    # which simplifies exactly to:
    #
    #   bd = b * (1-rflibc(minch))
    # ------------------------------------------------------------------

    depth = np.asarray(
        slope
        * reference_depth,
        dtype=np.float32,
    )

    # bandmp deletes both fit and depth if b is effectively zero.
    fit = np.ma.masked_where(
        invalid_slope,
        fit,
    )

    depth = np.ma.masked_where(
        invalid_slope,
        depth,
    )

    fit_depth = np.ma.asarray(
        fit * depth,
        dtype=np.float32,
    )

    # ------------------------------------------------------------------
    # Target continuum geometry at native minch/maxch.
    # ------------------------------------------------------------------

    continuum = np.ma.asarray(
        target.continuum,
        dtype=np.float32,
    )[..., extremum_index]

    left_continuum = np.ma.asarray(
        target.left_continuum,
        dtype=np.float32,
    )

    right_continuum = np.ma.asarray(
        target.right_continuum,
        dtype=np.float32,
    )

    band_bottom = np.ma.asarray(
        continuum
        * (
            np.float32(1.0)
            - depth
        ),
        dtype=np.float32,
    )

    reflectance_depth = np.ma.asarray(
        continuum * depth,
        dtype=np.float32,
    )

    # ------------------------------------------------------------------
    # Continuum/shoulder measurements.
    #
    # Kept in the existing form for now. Native denominator clamping
    # belongs to the feature-test implementation and is independent of
    # the current calcite/aragonite comparison.
    # ------------------------------------------------------------------

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        left_right_ratio = (
            left_continuum
            / right_continuum
        )

        right_left_ratio = (
            right_continuum
            / left_continuum
        )

        left_shoulder_ratio = (
            (
                left_continuum
                - band_bottom
            )
            /
            (
                right_continuum
                - band_bottom
            )
        )

        right_shoulder_ratio = (
            (
                right_continuum
                - band_bottom
            )
            /
            (
                left_continuum
                - band_bottom
            )
        )

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

        "value_count": reference_cr.count(),
    }


#==================== 

