"""
Pixel-independent preparation of every positive feature, and the per-strip feature fit built on it.

prepare_reference_features() does, once per run, everything in fit_feature that does not depend on pixels: the
status checks, the reference continuum removal (bdmset) and the reference half of characterise_feature. Entries are
keyed by reference_key(rule, feature), derived from the rule, so rules sharing a reference record and windows share
an entry. Each entry also carries the channel geometry of its target preparation (target_key); features with the
same target_key need the identical target continuum removal, so fit_target_group() fits them together:

    linear  one fused numba pass per pixel: bandmp window means, continuum, continuum-removed row, and the REAL*8
            least-squares sums for every reference in the group, in channel order as the native DO loops
    convex  NumPy: bandmpcv continuum (scipy CubicSpline, as curved_feature_continuum), then the sums batched

The results are what resolve_feature_tests reads from fit_feature: fit, depth, continuum at the band extremum,
left/right continuum means, polarity and reference_area. Arithmetic and thresholds are fit_feature's (REAL*4, with
float32 threshold comparisons); only the REAL*8 summation order can differ, at the 1 ulp level.
"""
from dataclasses import dataclass

import numpy as np
from numba import njit, prange
from scipy.interpolate import CubicSpline

from .config import REFERENCE_LIBRARIES
from .tetracorder_ops import (_CONT_MAX, _CONT_MIN, _TINY, _remove_continuum, _window_mean, _wtochbin,
                              linear_feature_continuum)


def reference_key(rule, feature):
    """Cache key of one positive feature: what its reference preparation depends on, derived from the rule."""
    record = rule["library_records"]["SMALL"]
    library = REFERENCE_LIBRARIES[record["library"].strip("[]")]
    windows = tuple(tuple(float(w) for w in window) for window in feature["windows"])
    return library, int(record["record"]), feature["continuum"], windows


@dataclass(frozen=True, eq=False)
class ReferenceFeature:
    status: str                                 # valid | out_of_range | disabled | invalid_reference | invalid_data
    continuum: str                              # linear | convex

    # geometry: pixel-independent, also defines the target preparation
    channels: tuple | None = None               # linear (cl1, cl2, cr1, cr2); convex n1..n8
    span: slice | None = None                   # fit and sum range: cl1..cr2 / n3..n6
    interior: np.ndarray | None = None          # bool over span, band channels only
    valid: np.ndarray | None = None             # bool over span: valid_bands & isfinite(reference)
    cwav: tuple | None = None                   # convex: window mid-wavelengths (spline anchors)
    target_key: tuple | None = None             # groups features whose target preparation is identical

    # reference side of bdmset / bandmp
    r: np.ndarray | None = None                 # float64 over span, continuum-removed reference, NaN -> 0
    rr: np.ndarray | None = None                # float64, r * r
    ref_ok: np.ndarray | None = None            # bool over span, isfinite(continuum-removed reference)
    extremum_index: int | None = None           # into span
    r_extremum: np.float32 | None = None        # rflibc[extremum_index]
    polarity: str | None = None
    reference_depth: np.float32 | None = None
    reference_area: float | None = None         # getifeat weight: sum |1 - cr| over interior

    


def _prepare_reference_feature(reference, feature, wavelengths, valid_bands):
    """
    One feature, in fit_feature's order so the status matches it: window coverage checks, reference continuum
    (invalid_reference), target geometry (invalid_data, convex only), then the reference checks of
    characterise_feature (invalid_reference).
    """
    continuum = feature["continuum"]
    windows = feature["windows"]
    valid = valid_bands & np.isfinite(np.asarray(reference, dtype=float))

    def failed(status):
        return ReferenceFeature(status=status, continuum=continuum)

    if continuum == "linear":
        (l0, l1), (r0, r1) = windows
        left = (wavelengths >= l0) & (wavelengths <= l1)
        right = (wavelengths >= r0) & (wavelengths <= r1)
        if not left.any() or not right.any():
            return failed("out_of_range")
        if not (left & valid).any() or not (right & valid).any():
            return failed("disabled")
        band = (wavelengths > l1) & (wavelengths < r0)
        reference_windows = windows[0], windows[1]
    elif continuum == "convex":
        for i, (start, stop) in enumerate(windows):
            selected = (wavelengths >= start) & (wavelengths <= stop)
            if not selected.any():
                return failed("out_of_range")
            if i in (1, 2) and not (selected & valid).any():
                return failed("disabled")
        band = (wavelengths > windows[1][1]) & (wavelengths < windows[2][0])
        reference_windows = windows[1], windows[2]          # bdmset: straight line through the inner windows
    else:
        raise ValueError(f"Unknown continuum type: {continuum}")

    if not band.any():
        return failed("out_of_range")
    if not (band & valid).any():
        return failed("disabled")

    try:
        ref = linear_feature_continuum(reference, wavelengths, *reference_windows, valid_bands=valid)
    except ValueError:
        return failed("invalid_reference")

    cwav = None
    if continuum == "convex":
        # the ValueError paths of curved_feature_continuum that depend only on geometry
        wl32 = np.asarray(wavelengths, dtype=np.float32)
        try:
            if len(windows) != 4:
                raise ValueError
            (n1, n2), (n3, n4), (n5, n6), (n7, n8) = (_wtochbin(wl32, float(a), float(b)) for a, b in windows)
            if n1 > n2 or n2 >= n3 or n3 > n4 or n4 + 1 >= n5 or n5 > n6 or n6 >= n7 or n7 > n8:
                raise ValueError
            mids = np.array([(float(a) + float(b)) / 2.0 for a, b in windows], dtype=np.float32)
            if not np.all(np.diff(mids) > 0):
                raise ValueError
        except ValueError:
            return failed("invalid_data")
        channels = (n1, n2, n3, n4, n5, n6, n7, n8)
        cwav = tuple(float(m) for m in mids)
    else:
        channels = tuple(int(c) for c in ref.channels)

    # reference half of characterise_feature
    rflibc = np.asarray(ref.continuum_removed, dtype=np.float32)
    interior = ref.interior
    inner = np.flatnonzero(interior & np.isfinite(rflibc))
    if inner.size == 0:
        return failed("invalid_reference")
    minch = inner[np.argmin(rflibc[inner])]
    maxch = inner[np.argmax(rflibc[inner])]
    rmin, rmax = rflibc[minch], rflibc[maxch]
    if rmin < 0.0:
        return failed("invalid_reference")
    if (rmax - 1.0) > (1.0 - rmin):
        polarity, extremum = "emission", int(maxch)
    else:
        polarity, extremum = "absorption", int(minch)

    ref_ok = np.isfinite(rflibc)
    r = np.where(ref_ok, rflibc, 0).astype(np.float64)
    cl1, _, _, cr2 = ref.channels
    span = slice(int(cl1), int(cr2) + 1)
    valid_span = valid[span]

    return ReferenceFeature(
        status="valid", continuum=continuum, channels=channels, span=span, interior=interior, valid=valid_span,
        cwav=cwav, target_key=(continuum, channels, cwav, valid_span.tobytes()),
        r=r, rr=r * r, ref_ok=ref_ok, extremum_index=extremum, r_extremum=rflibc[extremum], polarity=polarity,
        reference_depth=np.float32(1.0 - rflibc[extremum]),
        reference_area=float(np.nansum(np.abs(1.0 - rflibc[interior]))),
    )


def prepare_reference_features(rules, rule_spectra, wavelengths, valid_bands):
    """{reference_key: ReferenceFeature} for every positive feature of every tricorder-primary rule."""
    valid_bands = np.asarray(valid_bands, dtype=bool)
    cache = {}
    for rule_id, rule in rules.items():
        if rule["algorithm"] != "tricorder-primary":
            continue
        for feature in rule["features"]:
            if feature["role"] == "not":
                continue
            key = reference_key(rule, feature)
            if key not in cache:
                cache[key] = _prepare_reference_feature(rule_spectra[rule_id], feature, wavelengths, valid_bands)
    return cache


def group_by_target(reference_features):
    """{target_key: [reference_key, ...]}: references whose target preparation is identical."""
    groups = {}
    for key, entry in reference_features.items():
        if entry.status == "valid":
            groups.setdefault(entry.target_key, []).append(key)
    return groups


# ===================== convex groups: NumPy =====================

def _prepare_convex_target(spectra, wl32, entry):
    """
    Target side of curved_feature_continuum (bandmpcv) from the cached geometry: (continuum-removed (P, n),
    continuum (P, n), left mean (P,), right mean (P,)). Same arithmetic, same order.
    """
    n1, n2, n3, n4, n5, n6, n7, n8 = entry.channels
    means, counts, wl_means = [], [], []
    for i, (a, b) in enumerate(((n1, n2), (n3, n4), (n5, n6), (n7, n8))):
        vals = spectra[:, a:b + 1]
        if i in (1, 2):                     # valid bands honoured in the inner windows only (inside the span)
            vals = np.where(entry.valid[a - n3:b - n3 + 1], vals, np.nan).astype(np.float32)
        m, w, n = _window_mean(vals, wl32[a:b + 1], np.ones(b - a + 1, dtype=bool))
        means.append(m)
        counts.append(n)
        wl_means.append(w)
    anchors = np.stack(means, axis=-1).astype(np.float32)
    deleted = np.zeros(anchors.shape[:-1], dtype=bool)
    for m, n in zip(means, counts):
        deleted |= (n < 1) | ~(m >= _TINY)
    deleted |= ~(np.abs(wl_means[2] - wl_means[1]) >= _TINY)

    wav = wl32[entry.span]
    safe = np.where(deleted[:, None], np.float32(1.0), anchors)
    spline = CubicSpline(np.asarray(entry.cwav, dtype=float), safe.astype(float), axis=-1, bc_type="natural")
    continuum = np.asarray(spline(wav.astype(float)), dtype=np.float32)
    values = np.where(entry.valid, spectra[:, entry.span], np.nan).astype(np.float32)
    continuum_removed = _remove_continuum(values, continuum, np.ones(wav.shape, dtype=bool))
    continuum_removed[deleted] = np.nan
    return continuum_removed, continuum, means[1], means[2]


def _compare_batch(prep, entries):
    """
    characterise_feature for k references sharing one target preparation: (fit, depth, conref) as (P, k) and the
    left / right means as (P,). n, sumo and sumoo use characterise_feature's own NumPy calls (shared across the
    references with complete spans); the reference-dependent sums are (P, n) @ (n, k) products.
    """
    crem, cont, left_mean, right_mean = prep
    tok = np.isfinite(crem)
    tokf = tok.astype(np.float64)
    dro = np.where(tok, crem, 0).astype(np.float64)

    r = np.stack([e.r for e in entries], axis=1)
    rr = np.stack([e.rr for e in entries], axis=1)
    ref_ok = np.stack([e.ref_ok for e in entries], axis=1)

    suml = tokf @ r
    sumll = tokf @ rr
    sumol = dro @ r
    n = np.empty(suml.shape)
    sumo = np.empty(suml.shape)
    sumoo = np.empty(suml.shape)
    full = ref_ok.all(axis=0)
    if full.any():
        n[:, full] = tok.sum(axis=-1)[:, None]
        sumo[:, full] = dro.sum(axis=-1)[:, None]
        sumoo[:, full] = np.einsum("...i,...i->...", dro, dro)[:, None]
    for j in np.flatnonzero(~full):
        ok = tok & ref_ok[:, j]
        d = np.where(ok, crem, 0).astype(np.float64)
        n[:, j] = ok.sum(axis=-1)
        sumo[:, j] = d.sum(axis=-1)
        sumoo[:, j] = np.einsum("...i,...i->...", d, d)
    dxn = np.where(n > 0, n, 1).astype(np.float64)

    top = (sumol - sumo * suml / dxn).astype(np.float32)
    bottom = (sumll - suml * suml / dxn).astype(np.float32)
    botm2 = (sumoo - sumo * sumo / dxn).astype(np.float32)
    r_extremum = np.array([e.r_extremum for e in entries], dtype=np.float32)
    extremum = np.array([e.extremum_index for e in entries])

    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        slope = np.where(np.abs(bottom) < _TINY, 0, top / bottom).astype(np.float32)
        xk = ((1.0 - slope) / slope).astype(np.float32)
        xk1 = (xk + 1.0).astype(np.float32)
        rftemp = ((r_extremum + xk) / xk1).astype(np.float32)
        depth = (1.0 - rftemp).astype(np.float32)
        bprime = np.where(np.abs(botm2) < _TINY, 0, top / botm2).astype(np.float32)
        fit = np.sqrt(np.abs(slope * bprime)).astype(np.float32)

    conref = cont[:, extremum]
    deleted = ((n == 0) | ~(np.abs(slope) >= _TINY) | ~(np.abs(xk1) >= _TINY)
               | ~((fit > 0.0) & (fit < 1.1)) | ~(np.abs(depth) >= 0.1e-6))
    for c in (conref, left_mean[:, None], right_mean[:, None]):
        deleted |= ~((c >= _CONT_MIN) & (c <= _CONT_MAX))
    fit = np.where(deleted, np.nan, fit).astype(np.float32)
    depth = np.where(deleted, np.nan, depth).astype(np.float32)
    return fit, depth, conref, left_mean, right_mean


# ===================== linear groups: fused numba kernel =====================

# float32 thresholds: fit_feature compares float32 arrays with these in float32, so the kernel must too
_ONE32 = np.float32(1.0)
_ZERO32 = np.float32(0.0)
_NAN32 = np.float32(np.nan)
_FITMAX32 = np.float32(1.1)
_DEPTHMIN32 = np.float32(0.1e-6)
_CMIN32 = np.float32(_CONT_MIN)
_CMAX32 = np.float32(_CONT_MAX)

@njit(error_model="numpy", parallel=True, cache=True)
def _fused_linear(spectra, cl1, n_left, right_start, wav, valid, r, rr, ref_ok, r_extremum, extremum,
                  fit, depth, conref, left_out, right_out, crem):
    """bandmp + the characterise_feature sums for k references, one pixel at a time (see module docstring)."""
    n_pix = spectra.shape[0]
    nch = wav.shape[0]
    k = r.shape[0]
    for p in prange(n_pix):
        ls = _ZERO32
        lw = _ZERO32
        nl = 0
        for c in range(n_left):
            if valid[c]:
                v = spectra[p, cl1 + c]
                if np.isfinite(v):
                    ls += v
                    lw += wav[c]
                    nl += 1
        rs = _ZERO32
        rw = _ZERO32
        nr = 0
        for c in range(right_start, nch):
            if valid[c]:
                v = spectra[p, cl1 + c]
                if np.isfinite(v):
                    rs += v
                    rw += wav[c]
                    nr += 1
        lrefl = ls / np.float32(nl)                         # 0 / 0 -> NaN under error_model="numpy"
        lwl = lw / np.float32(nl)
        rrefl = rs / np.float32(nr)
        rwl = rw / np.float32(nr)
        bottom = rwl - lwl
        deleted = (nl < 1 or nr < 1 or not (lrefl >= _TINY) or not (rrefl >= _TINY)
                   or not (abs(bottom) >= _TINY))
        slope = (rrefl - lrefl) / bottom
        intercept = rrefl - slope * rwl
        left_out[p] = lrefl
        right_out[p] = rrefl

        for c in range(nch):
            cont = slope * wav[c] + intercept
            v = spectra[p, cl1 + c]
            if deleted or not valid[c] or not np.isfinite(v) or not (abs(cont) > _TINY):
                crem[p, c] = _NAN32
            else:
                crem[p, c] = v / cont

        for j in range(k):
            cnt = 0
            sl = 0.0
            sll = 0.0
            sol = 0.0
            so = 0.0
            soo = 0.0
            for c in range(nch):
                if ref_ok[j, c]:
                    o = crem[p, c]
                    if np.isfinite(o):
                        d = np.float64(o)
                        cnt += 1
                        sl += r[j, c]
                        sll += rr[j, c]
                        sol += d * r[j, c]
                        so += d
                        soo += d * d
            dxn = np.float64(cnt) if cnt > 0 else 1.0
            top = np.float32(sol - so * sl / dxn)
            bot = np.float32(sll - sl * sl / dxn)
            bot2 = np.float32(soo - so * so / dxn)
            if abs(bot) < _TINY:
                sj = _ZERO32
            else:
                sj = top / bot
            xk = (_ONE32 - sj) / sj
            xk1 = xk + _ONE32
            rftemp = (r_extremum[j] + xk) / xk1
            dep = _ONE32 - rftemp
            if abs(bot2) < _TINY:
                bprime = _ZERO32
            else:
                bprime = top / bot2
            ft = np.float32(np.sqrt(abs(sj * bprime)))
            cref = slope * wav[extremum[j]] + intercept
            bad = (cnt == 0 or not (abs(sj) >= _TINY) or not (abs(xk1) >= _TINY)
                   or not (ft > _ZERO32 and ft < _FITMAX32) or not (abs(dep) >= _DEPTHMIN32)
                   or not (cref >= _CMIN32 and cref <= _CMAX32)
                   or not (lrefl >= _CMIN32 and lrefl <= _CMAX32)
                   or not (rrefl >= _CMIN32 and rrefl <= _CMAX32))
            fit[p, j] = _NAN32 if bad else ft
            depth[p, j] = _NAN32 if bad else dep
            conref[p, j] = cref


def _fused_linear_group(spectra, wl32, entries, parallel):
    cl1, cl2, cr1, cr2 = entries[0].channels
    nch = cr2 - cl1 + 1
    n_pix, k = spectra.shape[0], len(entries)
    out = [np.empty((n_pix, k), dtype=np.float32) for _ in range(3)]                # fit, depth, conref
    left = np.empty(n_pix, dtype=np.float32)
    right = np.empty(n_pix, dtype=np.float32)
    _fused_linear(spectra, cl1, cl2 - cl1 + 1, cr1 - cl1,
           np.ascontiguousarray(wl32[cl1:cr2 + 1]), np.ascontiguousarray(entries[0].valid),
           np.ascontiguousarray(np.stack([e.r for e in entries])),
           np.ascontiguousarray(np.stack([e.rr for e in entries])),
           np.ascontiguousarray(np.stack([e.ref_ok for e in entries])),
           np.array([e.r_extremum for e in entries], dtype=np.float32),
           np.array([e.extremum_index for e in entries], dtype=np.int64),
           *out, left, right, np.empty((n_pix, nch), dtype=np.float32))
    return (*out, left, right)


# ===================== entry point =====================

def fit_target_group(spectra, wl32, entries, lead_shape=None, parallel=True):
    """
    Fit every reference of one target group. spectra: (pixels, bands) C-contiguous float32. Returns one dict per
    entry with the fit_feature keys resolve_feature_tests reads, arrays shaped lead_shape (default (pixels,)).
    """
    lead = (spectra.shape[0],) if lead_shape is None else tuple(lead_shape)
    if entries[0].continuum == "linear":
        fit, depth, conref, left, right = _fused_linear_group(spectra, wl32, entries, parallel)
    else:
        fit, depth, conref, left, right = _compare_batch(_prepare_convex_target(spectra, wl32, entries[0]), entries)
    left, right = left.reshape(lead), right.reshape(lead)
    return [{"fit": fit[:, j].reshape(lead), "depth": depth[:, j].reshape(lead),
             "continuum": conref[:, j].reshape(lead), "left_continuum": left, "right_continuum": right,
             "polarity": e.polarity, "reference_area": e.reference_area} for j, e in enumerate(entries)]