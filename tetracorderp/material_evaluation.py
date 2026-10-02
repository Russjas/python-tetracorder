"""
Class for evaluating a specific material from the tetracorder rule sets
"""
from pathlib import Path
import json
import numpy as np

from .tetracorder_ops import (linear_feature_continuum, 
                             curved_feature_continuum,
                             characterise_feature,
                             fuzzy_greater,
                             fuzzy_less,
                             continuum_test, 
                             load_references,
                             prepare_rules,
                             _wtochbin
)
from .config import RRATIO_REFERENCES
from .convolve import Convolver

def _zero_nan(x):
    """Deleted (NaN) pixels contribute 0, as np.ma.filled(x, 0.0) did for masked ones."""
    return np.where(np.isnan(x), 0.0, x)


class MaterialEvaluator:
    """
    Evaluate individual Tetracorder materials through:

        positive feature evaluation
            to
        feature role semantics
            to
        material fit / depth / fit-depth
            to
        material constraints
            to
        NOT vetoes

    Group and case resolution are intentionally outside this class.
    """

    def __init__(self,
                 rules: dict|str|Path,
                 references: dict|str|Path,
                 target_wavelengths: np.ndarray,
                 target_fwhm: np.ndarray,
                 target_valid_bands: np.ndarray|None = None,
                 mode: str = "default",
                 disabled_materials: set|None = None
                 ):

        self.rules, self.not_definitions = self._prepare_rules(rules, mode)
        self.disabled_materials = set(disabled_materials or ())
        self._not_source_cache = {}

        self.wavelengths = target_wavelengths
        self.fwhm = target_fwhm
        self.valid_bands = target_valid_bands

        if self.valid_bands is None: 
            self.valid_bands = np.ones(self.wavelengths.shape, dtype=bool)
        else:
            self.valid_bands = np.asarray(self.valid_bands, dtype=bool)
        if self.valid_bands.shape != self.wavelengths.shape:
            raise ValueError("valid_bands must have the same shape as target_wavelengths")
        if isinstance(references, (str, Path)):
            references = load_references(references)
        self.rule_spectra = self._prepare_rule_spectra(references)
        self.rratio_spectra = self._prepare_rratio_spectra(references)


    def cache_clear(self):
        """Forget cached NOT-source results. Call before evaluating a new set of spectra with same_target=True."""
        self._not_source_cache.clear()


    def evaluate(self, rule_id: str, target_spectra: np.ndarray, same_target: bool = False) -> dict | None:

        if not same_target:
            self.cache_clear()

        if rule_id in self.disabled_materials:
            return None
        
        rule = self.rules[rule_id]

        if rule["algorithm"] == "nvres":
            return self._evaluate_nvres(rule_id, rule, target_spectra)
        
        if rule["algorithm"] != "tricorder-primary":
            return None

        rule_features = {feature["id"]: feature for feature in rule["features"]}

        reference_spectrum = self.rule_spectra[rule_id]

        resolved_features = self._resolve_positive_features(rule_features, reference_spectrum,
                                                            target_spectra)

        material = resolve_positive_material(rule_features, resolved_features)

        if material is None: 
            return None

        material = resolve_material_constraints(material, rule.get("constraints", []))

        material = self._resolve_nots(material, rule_features, resolved_features, target_spectra)

        return material

    def _evaluate_nvres(self, rule_id, rule, target_spectra,):
        """
        Evaluate one of the nvres case rules.
        """

        if rule["kind"] != "case":
            raise ValueError(f"Unexpected non-case nvres rule: {rule_id}")

        if len(rule["features"]) != 1:
            raise ValueError(f"nvres rule {rule_id} does not have exactly one feature")

        feature = rule["features"][0]

        # ----------------------------------------------------------
        # Resolve RRATIO support spectrum.
        # ----------------------------------------------------------

        preratio = rule.get("preratio")

        if preratio is None:
            raise ValueError(f"nvres rule {rule_id} has no RRATIO")

        ratio_type, ratio_name = preratio.split()

        if ratio_type != "RRATIO":
            raise ValueError(f"Unexpected nvres preratio: {preratio}")

        if ratio_name not in self.rratio_spectra:
            raise ValueError(f"reference spectra for {ratio_name} was not initialised")       

        ratio_spectrum = self.rratio_spectra[ratio_name]

        # Main 7632 / 7650 reference is already prepared normally.
        reference_spectrum = self.rule_spectra[rule_id]

        feature_result = resolve_nvres_feature(
            reference_spectrum,
            ratio_spectrum,
            feature,
            target_spectra,
            self.wavelengths,
            self.valid_bands,
        )

        rule_features = {feature["id"]: feature}

        resolved_features = {feature["id"]: feature_result}

        material = resolve_positive_material(rule_features, resolved_features)

        if material is None:
            return None

        material = resolve_material_constraints( material, rule.get("constraints", []))

        return material


    def _prepare_rules(self, rules, mode):
        if isinstance(rules, (str, Path)):
            with open(rules, "r", encoding="utf-8") as f:
                rules = json.load(f)

        rules = prepare_rules(rules, mode=mode)

        return rules["rules"], rules["not_definitions"]


    def _prepare_rule_spectra(self, references):
        LIBRARY_MAP = {
            "splib06": "splib06b",
            "sprlb06": "sprlb06b",
        }

        if isinstance(references, (str, Path)):
            references = load_references(references)

        spectra = {}
        cache = {}
        fwhm = np.broadcast_to(np.asarray(self.fwhm, dtype=float), self.wavelengths.shape)
        convolver = Convolver(self.wavelengths, fwhm)

        for rule_id, rule in self.rules.items():
            library = rule["library_records"]["SMALL"]["library"].strip("[]")
            library = LIBRARY_MAP[library]

            record = int(rule["library_records"]["SMALL"]["record"])
            key = (library, record)

            if key not in cache:
                reference = references[key]
                cache[key] = convolver.convolve(reference["wavelengths"],
                                                reference["fwhm"],
                                                reference["reflectance"])
            spectra[rule_id] = cache[key]
        
        return spectra

    def _prepare_rratio_spectra(self, references):
        spectra = {}
        for ratio_name, key in RRATIO_REFERENCES.items():
            ratio_reference = references[key]
            fwhm = np.broadcast_to(np.asarray(self.fwhm, dtype=float), self.wavelengths.shape)
            convolver = Convolver(self.wavelengths, fwhm)

            spectra[ratio_name] = convolver.convolve(ratio_reference["wavelengths"],
                                                    ratio_reference["fwhm"],
                                                    ratio_reference["reflectance"])
        return spectra



    def _resolve_positive_features(self,
                                   rule_features: dict,
                                   reference_spectrum: np.ndarray,
                                   target_spectra: np.ndarray) -> dict:

        resolved_features = {}

        for feature_id, feature in rule_features.items():
            if feature["role"] == "not":
                continue

            resolved_features[feature_id] = resolve_feature_tests(reference_spectrum,
                                            feature, target_spectra,
                                            self.wavelengths, self.valid_bands)

        return resolved_features


    def _native_not_source(self, source_rule_id, source_feature_id, target_spectra):
        """
        Fit and depth of a NOT source feature as native tp1mat leaves them in
        zfit / zdepth. Cached per (source_rule_id, source_feature_id): many
        materials in a ruleset veto against the same common source.

        The result also depends on target_spectra, which is deliberately not in
        the key. Cache validity is a contract with the caller: evaluate() clears
        the cache unless called with same_target=True, a promise that the pixels
        are unchanged since the last cache_clear().
        """
        key = (source_rule_id, source_feature_id)
        if key in self._not_source_cache:
            return self._not_source_cache[key]

        result = self._native_not_source_uncached(source_rule_id, source_feature_id, target_spectra)
        self._not_source_cache[key] = result
        return result

    def _native_not_source_uncached(self, source_rule_id, source_feature_id, target_spectra):
        """
        tp1mat evaluates a material's features in rule order. When a weak,
        diagnostic or must-have feature is absent (depth / sign <= 1e-6) it
        zeroes that feature and every later feature, then returns. So the
        source feature is 0 wherever an earlier W/D/M feature of the source
        material failed. Disabled features (ifeatenable = 0) do not trigger
        this.
        """
        source_rule = self.rules[source_rule_id]
        reference = self.rule_spectra[source_rule_id]
        dead = None

        positive = [f for f in source_rule["features"] if f["role"] != "not"]
        # native tp1mat loops ifeat = 1..nfeat, ifeat taken from the label (f2a -> 2)
        for feature in sorted(positive, key=lambda f: f["number"]):
            if feature["role"] == "not":
                continue

            result = resolve_feature_tests(reference, feature, target_spectra,
                                           self.wavelengths, self.valid_bands)
            valid = result["status"] == "valid"

            if valid:
                fit = _zero_nan(np.asarray(result["mod_fit"], dtype=float))
                depth = _zero_nan(np.asarray(result["mod_depth"], dtype=float))
                if dead is None:
                    dead = np.zeros(fit.shape, dtype=bool)
                if feature["role"] in {"weak", "diagnostic", "must_have"}:
                    sign = 1.0 if result["polarity"] == "absorption" else -1.0
                    with np.errstate(invalid="ignore"):
                        dead_here = dead | ~(depth / sign > 1e-6)
                else:
                    dead_here = dead

            if feature["id"] == source_feature_id:
                if not valid:
                    return result
                return {**result,
                        "mod_fit": np.where(dead_here, 0.0, fit),
                        "mod_depth": np.where(dead_here, 0.0, depth)}

            if valid:
                dead = dead_here

        raise KeyError(f"{source_rule_id} has no feature {source_feature_id}")


    def _resolve_nots(self,
                      material_result: dict,
                      rule_features: dict,
                      resolved_features: dict,
                      target_spectra: np.ndarray) -> dict:

        fit = material_result["fit"].copy()
        depth = material_result["depth"].copy()
        fit_depth = material_result["fit_depth"].copy()

        alive = fit > 0

        if not np.any(alive):
            return {
                "fit": fit,
                "depth": depth,
                "fit_depth": fit_depth,
                "rejected": material_result["rejected"],
            }

        for feature_id, not_feature in rule_features.items():
            if not_feature["role"] != "not":
                continue

            # A disabled material is never evaluated, so native's stored
            # zdepth/zfit for it stay zero and the NOT silently never fires.
            if not_feature["source_material"] in self.disabled_materials:
                continue
            
            source_result = self._native_not_source(not_feature["source_material"],
                                                    not_feature["source_feature"],
                                                    target_spectra)

            not_mask = resolve_not(not_feature, source_result, resolved_features)

            if not_mask is None:
                continue

            reject = alive & not_mask

            fit[reject] = 0.0
            depth[reject] = 0.0
            fit_depth[reject] = 0.0

            alive[reject] = False

            if not np.any(alive):
                break

        return {
            "fit": fit,
            "depth": depth,
            "fit_depth": fit_depth,
            "rejected": material_result["rejected"],
        }


    def evaluate_all(self, target_spectra: np.ndarray) -> dict:
        results = {}
        self.cache_clear()
        for rule_id in self.rules:
            results[rule_id] = self.evaluate(rule_id, target_spectra, same_target=True)

        return results



def fit_feature(reference, target, target_wavelengths, windows, continuum = "linear", 
                valid_bands = None, polarity = "absorption"):
    """
    Fits a feature using python equivalents of the specpr operations

    Parameters
    ----------
    reference : ndarray
        reference reflectance spectrum convolved with target wavelengths
    target : ndarray
        Target spectrum or cube
    target_wavelengths : ndarray
        Wavelengths of the target spectrum or cube.
    windows : list
         Feature start and end windows to define the continuum, specified in the rule 
    continuum : str, optional
        "linear" or "convex". Continuum calculation type specified in the rule
    valid_bands : numpy.ndarray of bool, optional
        Usable bands. Bands where the reference is deleted (NaN) are dropped as well.
    polarity : {"absorption", "emission"}, default "absorption"
        Passed to characterise_feature, which takes the feature type from the reference regardless.


    Returns
    -------
    result : dict or None
        characterise_feature output, or None when the feature cannot be fitted.
    status : str
        "valid", or why the feature was not fitted: "out_of_range" (a window outside the sensor), "disabled"
        (no usable band in a window), "invalid_reference" or "invalid_data".

    Raises
    ------
    ValueError
        on unknown continuum type.   
    """
    if valid_bands is None:
        valid_bands = np.ones(target_wavelengths.shape, dtype=bool)
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)
    if valid_bands.shape != target_wavelengths.shape:
        raise ValueError("valid_bands must have the same shape as target_wavelengths")

    # bandmp skips a channel wherever the reference is deleted (rflibc == delpt)
    valid_bands = valid_bands & np.isfinite(np.asarray(reference, dtype=float))
    
    if continuum == "linear":
        left_window = windows[0]
        right_window = windows[1]
        
        left_selected = ((target_wavelengths >= left_window[0])
            & (target_wavelengths <= left_window[1]))

        right_selected = ((target_wavelengths >= right_window[0])
            & (target_wavelengths <= right_window[1]))

        if not np.any(left_selected) or not np.any(right_selected):
            return None, "out_of_range"

        # Check that each continuum window contains at least one usable band.
        if (not np.any(left_selected & valid_bands) or not np.any(right_selected & valid_bands)):
            return None, "disabled"

        # Check the actual feature interval separately.
        feature_selected = ((target_wavelengths > left_window[1]) & (target_wavelengths < right_window[0]))
        if not np.any(feature_selected):
            return None, "out_of_range"

        if not np.any(feature_selected & valid_bands):
            return None, "disabled"
       
        try:
            ref_cont = linear_feature_continuum(reference, target_wavelengths,
                                            left_window,right_window,valid_bands=valid_bands)
        except ValueError:
            return None, "invalid_reference"
        try:
            target_cont = linear_feature_continuum(target, target_wavelengths,
                                               left_window, right_window,valid_bands=valid_bands)
        except ValueError:
            return None, "invalid_data"
        
    elif continuum == "convex":
       
        for i, (start, stop) in enumerate(windows):
            selected = ((target_wavelengths >= start) & (target_wavelengths <= stop))

            if not np.any(selected):
                return None, "out_of_range"

            # native bandmpcv skips reference-deleted channels only in the inner
            # windows (rlbc is 0.0 outside them), so only those need one.
            if i in (1, 2) and not np.any(selected & valid_bands):
                return None, "disabled"
        
        feature_selected = ((target_wavelengths > windows[1][1])
                            & (target_wavelengths < windows[2][0]))

        if not np.any(feature_selected):
            return None, "out_of_range"

        if not np.any(feature_selected & valid_bands):
            return None, "disabled"
        
        try:
            # Tetracorder/Specpr: bdmset linear continuum on the inner windows for the reference
            # deliberate asymmetry between reference and target continuum
            ref_cont = linear_feature_continuum(reference, target_wavelengths,
                                                windows[1], windows[2],
                                                valid_bands=valid_bands)
        except ValueError:
            return None, "invalid_reference"
        try:
            target_cont = curved_feature_continuum(target, target_wavelengths,
                                               windows,valid_bands=valid_bands)
        except ValueError:
            return None, "invalid_data"

    else:
        raise ValueError(f"Unknown continuum type: {continuum}")
    try:
        result = characterise_feature(ref_cont, target_cont, polarity = polarity)
    except ValueError:
        return None, "invalid_reference"

    # reference_area (feature weight) comes from characterise_feature:
    # native getifeat sums |1 - cr| over the band channels only
    return result, "valid"


def fit_nvres_feature(
    reference,
    ratio_reference,
    target,
    target_wavelengths,
    windows,
    valid_bands=None,
):
    """
    Python implementation of the Tetracorder nvres red-edge
    feature calculation.

    reference
        Main rule reference spectrum, convolved to target wavelengths.

    ratio_reference
        RRATIO reference spectrum, convolved to target wavelengths.

    target
        Target spectrum or cube.

    Returns
    -------
    result, status

    result has the same basic structure as fit_feature(), but
    result["depth"] is the nvres normalized depth (bdnorm).
    """

    if valid_bands is None:
        valid_bands = np.ones(target_wavelengths.shape, dtype=bool)
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)

    left_window = windows[0]
    right_window = windows[1]

    left_selected = ((target_wavelengths >= left_window[0])
        & (target_wavelengths <= left_window[1]))

    right_selected = ((target_wavelengths >= right_window[0])
        & (target_wavelengths <= right_window[1]))

    feature_selected = ((target_wavelengths > left_window[1])
        & (target_wavelengths < right_window[0]))

    # Distinguish sensor coverage from deliberately disabled bands.
    if (not np.any(left_selected)
        or not np.any(right_selected)
        or not np.any(feature_selected)):
        return None, "out_of_range"

    if (not np.any(left_selected & valid_bands)
        or not np.any(right_selected & valid_bands)
        or not np.any(feature_selected & valid_bands)):
        return None, "disabled"

    # specpr wtochbin takes a CONTIGUOUS channel run and stops at the first
    # channel past w2. The AVIRIS array folds back at the spectrometer joins
    # (ch 32 -> 33 is 0.687 -> 0.664), so a wavelength mask picks up channels
    # from the far side that native never sees: 0.661-0.681 is channels 30-31
    # natively, 30, 31, 33, 34 by mask.
    cl1, cl2 = _wtochbin(target_wavelengths, left_window[0], left_window[1])
    cr1, cr2 = _wtochbin(target_wavelengths, right_window[0], right_window[1])

    region = slice(cl1, cr2 + 1)

    wavelengths = target_wavelengths[region]
    reference = reference[region]
    ratio_reference = ratio_reference[region]
    target = np.asarray(target)[..., region]
    valid = valid_bands[region]

    ratio_valid = (valid
        & np.isfinite(ratio_reference)
        & (np.abs(ratio_reference) > 1e-13))

    left = np.zeros(wavelengths.shape, dtype=bool)
    left[: cl2 - cl1 + 1] = True
    left &= ratio_valid

    right = np.zeros(wavelengths.shape, dtype=bool)
    right[cr1 - cl1: cr2 - cl1 + 1] = True
    right &= ratio_valid

    if not np.any(left) or not np.any(right):
        return None, "invalid_reference"

    # ----------------------------------------------------------
    # Continuum-window means.
    #
    # Ratfor:
    #   avlcv = left mean, RRATIO reference
    #   avlcu = left mean, unknown/target
    #   avrcv = right mean, RRATIO reference
    #   avrcu = right mean, unknown/target
    #
    # The RRATIO mean must use the same surviving target channels.
    # ----------------------------------------------------------

    def paired_means(selected):
        
        target_window = target[..., selected]
        ratio_window = np.broadcast_to(ratio_reference[selected], target_window.shape)
        ok = np.isfinite(target_window) & np.isfinite(ratio_window)
        count = ok.sum(axis=-1)
        # MaskedArray.mean: float32 sum of filled(0), then * 1. / count -> float64; empty -> NaN
        with np.errstate(divide="ignore", invalid="ignore"):
            target_mean = np.where(ok, target_window, np.float32(0)).sum(axis=-1) * 1. / count
            ratio_mean = np.where(ok, ratio_window, np.float32(0)).sum(axis=-1) * 1. / count
        return target_mean, ratio_mean

    avlcu, avlcv = paired_means(left)
    avrcu, avrcv = paired_means(right)

    # ----------------------------------------------------------
    # Normalize the RRATIO spectrum to the observed red edge.
    #
    # Ratfor:
    #
    #   a = (avlcu*avrcv - avrcu*avlcv)
    #       / (avrcu - avlcu)
    #
    #   b = (avlcv + a) / avlcu
    #
    #   vegspec = (vegspec + a) / b
    # ----------------------------------------------------------

    with np.errstate(divide="ignore", invalid="ignore"):

        a = ((avlcu * avrcv) - (avrcu * avlcv)) / (avrcu - avlcu)

        b = (avlcv + a) / avlcu

        normalized_ratio = (ratio_reference + a[..., None]) / b[..., None]

        ratioed_target = (target / normalized_ratio)

    ratioed_target = np.where(np.isfinite(ratioed_target), ratioed_target, np.nan)

    # RRATIO-deleted channels are also unusable by bandmp.
    ratioed_target[..., ~ratio_valid] = np.nan

    # ----------------------------------------------------------
    # Ordinary Specpr band matching on the ratioed spectrum.
    #
    # These nvres rules are emission features. Ratfor nftype = -1.
    # ----------------------------------------------------------

    result, status = fit_feature(
        reference,
        ratioed_target,
        wavelengths,
        windows,
        continuum="linear",
        valid_bands=ratio_valid,
        polarity="emission",
    )

    if status != "valid":
        return None, status

    # Preserve ordinary bandmp depth for debugging/parity testing.
    band_depth = result["depth"].copy()

    # ----------------------------------------------------------
    # Normalize the band depth for red-edge strength.
    #
    # Ratfor:
    #
    #   xndvivegspec = (avrcv-avlcv)/(avrcv+avlcv)
    #   xndviobs     = (avrcu-avlcu)/(avrcu+avlcu)
    #   xfactor      = xndvivegspec/xndviobs
    #   bdnorm       = bdepth * xfactor
    # ----------------------------------------------------------

    with np.errstate(divide="ignore", invalid="ignore"):

        reference_edge = ((avrcv - avlcv) / (avrcv + avlcv))
        target_edge = ((avrcu - avlcu) / (avrcu + avlcu))

        depth_factor = (reference_edge / target_edge)

        normalized_depth = (band_depth * depth_factor)

    normalized_depth = np.where(np.isfinite(normalized_depth), normalized_depth, np.nan)

    result["nvres_band_depth"] = band_depth
    result["nvres_depth_factor"] = depth_factor

    # nvres1mat uses bdnorm, not ordinary bdepth.
    result["depth"] = normalized_depth
    result["fit_depth"] = (result["fit"] * normalized_depth)

    return result, "valid"


def resolve_nvres_feature(
    reference_spectrum,
    ratio_spectrum,
    feature,
    target_spectra,
    target_wavelengths,
    valid_bands,
):
    """
    nvres1mat feature-level processing for the nvres rules
    present in the current Tetracorder rule set.
    """
    result, status = fit_nvres_feature(
        reference_spectrum,
        ratio_spectrum,
        target_spectra,
        target_wavelengths,
        feature["windows"],
        valid_bands=valid_bands,
    )

    if status != "valid":
        return {
            "status": status,
            "mod_fit": None,
            "mod_depth": None,
            "mod_fit_depth": None,
            "raw_fit": None,
            "raw_depth": None,
            "polarity": None,
            "reference_area": None,
        }

    fit = np.asarray(result["fit"]).copy()

    depth = np.asarray(result["depth"]).copy()

    # Current nvres rules contain only ct.
    for test in feature.get("tests", []):

        if test["test"] != "ct":
            raise NotImplementedError(f"Unexpected nvres feature test: {test['test']}")

        values = test["values"]

        if (isinstance(values, (list, tuple)) and len(values) == 1):
            values = values[0]

        passed = continuum_test(result["continuum"], values)

        fail = ~passed

        fit[fail] = 0.0
        depth[fail] = 0.0

    return {
        "status": "valid",
        "mod_fit": fit,
        "mod_depth": depth,
        "mod_fit_depth": fit * depth,

        # In nvres terms this is already bdnorm.
        "raw_fit": result["fit"],
        "raw_depth": result["depth"],

        "polarity": result["polarity"],
        "reference_area": result["reference_area"],
    }


def resolve_feature_tests(reference_spectrum, feature_dict, target_spectra, target_wvls, valid_bands = None):
    """
    Fit one feature and apply its tests in tp1mat order: the hard continuum tests (ct, lct, rct) zero failing
    pixels, then the fuzzy tests (lct/rct>, rct/lct>, rcbblc, lcbbrc, r*bd>) scale fit and depth, each applied
    to the current values.

    Parameters
    ----------
    reference_spectrum : numpy.ndarray
        Rule reference convolved to the target bands, shape (bands,).
    feature_dict : dict
        One feature of a prepared rule: windows, continuum and tests.
    target_spectra : numpy.ndarray
        Target spectra, (..., bands).
    target_wvls : numpy.ndarray
        Band centres in µm, shape (bands,).
    valid_bands : numpy.ndarray of bool, optional
        Usable bands.

    Returns
    -------
    dict
        status (as fit_feature); mod_fit, mod_depth, mod_fit_depth after the tests; raw_fit, raw_depth before
        them; polarity; reference_area (the feature weight). The arrays are None unless status is "valid".
    """
    
    windows = feature_dict["windows"]
    continuum_type = feature_dict["continuum"]
    feat_fit, status = fit_feature(reference_spectrum, target_spectra, target_wvls, windows, 
                                   continuum = continuum_type, valid_bands = valid_bands)
    if status != "valid":
        return {
            "status": status,
            "mod_fit": None,
            "mod_depth": None,
            "mod_fit_depth": None,
            "raw_fit": None,
            "raw_depth": None,
            "polarity": None,
            "reference_area": None,
        }

    
    
    
    tests = feature_dict.get("tests", [])

    fit = np.asarray(feat_fit["fit"]).copy()
    depth = np.asarray(feat_fit["depth"]).copy()

    continuum = np.asarray(feat_fit["continuum"])
    left_continuum = np.asarray(feat_fit["left_continuum"])
    right_continuum = np.asarray(feat_fit["right_continuum"])

    # Convenient lookup. A feature should only have one of each test type.
    tests_by_name = {test["test"]: test["values"] for test in tests}

    test_results = {}

    def values_for(name):
        values = tests_by_name[name]

        if isinstance(values, (list, tuple)) and len(values) == 1:
            return values[0]

        return values

    # ----------------------------------------------------------
    # 1. Hard continuum tests
    # Here any fails have fit and depth immediately zerod
    # ----------------------------------------------------------

    if "ct" in tests_by_name:
        result = continuum_test(
            continuum, values_for("ct"))
        test_results["ct"] = result
        fail = ~result
        fit[fail] = 0.0
        depth[fail] = 0.0

    if "lct" in tests_by_name:
        result = continuum_test(left_continuum, values_for("lct"))
        test_results["lct"] = result

        fail = ~result
        fit[fail] = 0.0
        depth[fail] = 0.0

    if "rct" in tests_by_name:
        result = continuum_test(right_continuum, values_for("rct"))
        test_results["rct"] = result

        fail = ~result
        fit[fail] = 0.0
        depth[fail] = 0.0

    # ----------------------------------------------------------
    # Helper for fuzzy attenuation.
    #
    # Tetracorder applies each fuzzy result immediately to the
    # CURRENT fit and depth.
    # ----------------------------------------------------------

    def apply_factor(name, factor):
        test_results[name] = factor

        fit[...] = fit * factor
        depth[...] = depth * factor

    # ----------------------------------------------------------
    # 2. Left continuum / right continuum
    # ---------------------------------------------------------

    if "lct/rct>" in tests_by_name:
        with np.errstate(divide="ignore", invalid="ignore"):
            value = left_continuum / right_continuum

        factor = fuzzy_greater(value, values_for("lct/rct>"))

        apply_factor("lct/rct>", factor)

    # ---------------------------------------------------------
    # 3. Right continuum / left continuum
    # ----------------------------------------------------------

    if "rct/lct>" in tests_by_name:
        with np.errstate(divide="ignore", invalid="ignore"):
            value = right_continuum / left_continuum

        factor = fuzzy_greater(value, values_for("rct/lct>"))

        apply_factor("rct/lct>", factor)

    # ----------------------------------------------------------
    # From this point onwards, band bottom must be calculated
    # from the CURRENT depth, not characterise_feature()["band_bottom"].
    # ----------------------------------------------------------

    def current_band_bottom():
        return continuum * (1.0 - depth)

    def current_right_shoulder_ratio():
        band_bottom = current_band_bottom()

        with np.errstate(divide="ignore", invalid="ignore"):
            return ((right_continuum - band_bottom) / (left_continuum - band_bottom))

    def current_left_shoulder_ratio():
        band_bottom = current_band_bottom()

        with np.errstate(divide="ignore", invalid="ignore"):
            return ((left_continuum - band_bottom) / (right_continuum - band_bottom))

    # ----------------------------------------------------------
    # 4. (right continuum - band bottom) /
    #    (left continuum - band bottom)
    # ----------------------------------------------------------

    if "rcbblc>" in tests_by_name:
        value = current_right_shoulder_ratio()

        factor = fuzzy_greater(value, values_for("rcbblc>"))

        apply_factor("rcbblc>", factor)

    # Recalculate because depth may just have changed.
    if "rcbblc<" in tests_by_name:
        value = current_right_shoulder_ratio()

        factor = fuzzy_less(value, values_for("rcbblc<"))

        apply_factor("rcbblc<", factor)

    # ----------------------------------------------------------
    # 5. (left continuum - band bottom) /
    #    (right continuum - band bottom)
    # ----------------------------------------------------------

    if "lcbbrc>" in tests_by_name:
        value = current_left_shoulder_ratio()

        factor = fuzzy_greater(value,values_for("lcbbrc>"))

        apply_factor("lcbbrc>", factor)

    # Again recalculate from the newly modified depth.
    if "lcbbrc<" in tests_by_name:
        value = current_left_shoulder_ratio()

        factor = fuzzy_less(value, values_for("lcbbrc<"))

        apply_factor("lcbbrc<", factor)

    # ----------------------------------------------------------
    # 6. Reflectance * |band depth|
    #
    # This also uses CURRENT depth after all previous attenuation.
    # ----------------------------------------------------------

    
    # tp1mat.r takes the r*bd> reject (lower) limit from zcontlgtr(1) (getifeat.r:506/851), i.e. the lct/rct>
    # lower limit, not from zrtimesbd(1), which is parsed (getifeat.r:1185) but never read:
    #   - no lct/rct> on the feature -> reject limit 0.0 (811 features)
    #   - lct/rct> present           -> reject limit = that test's lower ratio limit (10 features)
    # Reproduced here for parity with native: thresholds have been tuned around this behaviour, so changing it
    # is for upstream to decide (spectroscopy-tetracorder issue #6). The commented-out block below uses the
    # rule's own r*bd> limits instead.
#================un-comment this block below to end bug emulation =================================================    
    
    #if "r*bd>" in tests_by_name:
        #sign = 1.0 if feat_fit["polarity"] == "absorption" else -1.0
        #with np.errstate(invalid="ignore"):
        #    present = (depth / sign) > 1e-6
        #value = np.where(present, continuum * depth, 0.0)

        #factor = fuzzy_greater(value, values_for("r*bd>"))

        #apply_factor("r*bd>", factor)
#=================================================================================================================
#================comment the block below to end bug emulation =====================================================
    if "r*bd>" in tests_by_name:
        # tp1mat.r: xtmp1 = conref*bdepth - SIGNED, despite the Ratfor comment
        # claiming r*abs(bd) - and forced to 0.0 when the feature is absent
        # (bdepth/xfeat <= 0.1e-5). An emission feature therefore always gives
        # xtmp1 <= 0 and can never clear a lower limit of 0.0, so native never
        # detects a material whose emission feature carries an r*bd> test.
        sign = 1.0 if feat_fit["polarity"] == "absorption" else -1.0
        with np.errstate(invalid="ignore"):
            present = (depth / sign) > 1e-6
        value = np.where(present, continuum * depth, 0.0)

        # native limits (see above): reject from the lct/rct> lower limit, else 0.0; full from the rule
        rbd_upper = values_for("r*bd>")[1]
        #rbd_lower = values_for("r*bd>")[0]
        if "lct/rct>" in tests_by_name:
            rbd_lower = values_for("lct/rct>")[0]
        else:
            rbd_lower = 0.0

        factor = np.ones(np.shape(value), dtype=float)

        reject = value < rbd_lower
        fuzzy = (
            (value >= rbd_lower)
            & (value < rbd_upper)
        )

        factor[reject] = 0.0

        factor[fuzzy] = (
            (value[fuzzy] - rbd_lower)
            / (rbd_upper - rbd_lower)
        )

        apply_factor("r*bd>", factor)
#==========================================================================================
    # ----------------------------------------------------------
    # "weight" is not a spectral rejection/fuzzy test here.
    # Deal with it during feature weighting/material aggregation.
    # ----------------------------------------------------------

    if "weight" in tests_by_name:
        test_results["weight"] = values_for("weight")

    return {
        "status": "valid",
        "mod_fit": fit,
        "mod_depth": depth,
        "mod_fit_depth": fit * depth,
        "raw_fit": feat_fit["fit"],
        "raw_depth": feat_fit["depth"],
        "polarity": feat_fit["polarity"],
        "reference_area": feat_fit["reference_area"],
    }



def resolve_feature_weights(rule_features, resolved_features):
    raw_weights = {}

    for f_id, feature in rule_features.items():
        if feature["role"] == "not":
            continue

        result = resolved_features[f_id]

        if feature["role"] == "weak" or result["status"] != "valid":
            raw_weights[f_id] = 0.0
            continue

        weight_modifier = 1.0

        for test in feature.get("tests", []):
            if test["test"] == "weight":
                values = test["values"]
                weight_modifier = values[0] if isinstance(values, (list, tuple)) else values
                break

        raw_weights[f_id] = result["reference_area"] * weight_modifier

    total = sum(raw_weights.values())

    if total <= 0.0:
        return {f_id: 0.0 for f_id in raw_weights}

    return {f_id: weight / total for f_id, weight in raw_weights.items()}


def resolve_positive_material(rule_features, resolved_features):
    """
    Resolve material fit, depth and fit-depth from positive features.
    
    Role semantics - Preserved from Tetracorder documentation
    O = optional
    W = weak must-be-present
    D = diagnostic
    M = must-have diagnostic

    role         out_of_range / disabled        valid but fit/depth fails        valid and survives feature tests

    optional     ignore                          contributes zero                 contributes
    diagnostic   ignore                          reject material                  contributes
    weak         reject material                 reject material                  DOES NOT contribute to weighted sums
    must_have    reject material                 reject material                  contributes

    """ 
    

    first_valid = next((resolved_features[f_id] for f_id, feature in rule_features.items()
            if feature["role"] != "not" and resolved_features[f_id]["status"] == "valid"),
            None)
    if first_valid is None:
        return None

    feature_weights = resolve_feature_weights(rule_features, resolved_features)
    
    shape = np.shape(first_valid["mod_fit"])

    sum_fit = np.zeros(shape, dtype=float)
    sum_depth = np.zeros(shape, dtype=float)
    sum_fit_depth = np.zeros(shape, dtype=float)

    rejected = np.zeros(shape, dtype=bool)

    for f_id, feature in rule_features.items():
        if feature["role"] == "not":
            continue

        result = resolved_features[f_id]
        role = feature["role"]

        # Ratfor ifeatenable == 0.
        if result["status"] != "valid":
            if role in {"weak", "must_have"}:
                rejected[...] = True

            continue

        fit = np.asarray(result["mod_fit"])
        depth = np.asarray(result["mod_depth"])

        feature_sign = 1.0 if result["polarity"] == "absorption" else -1.0
        weight = feature_weights[f_id]

        # Tetracorder Ratfor:
        # xx = 1.0
        # if (bdepth / xfeat <= 0.1e-5) xx = 0.0
        with np.errstate(invalid="ignore"):
            present = (depth / feature_sign) > 1e-6        # deleted (NaN) -> False
        xx = present.astype(float)

        # featimprt > 0 means W, D or M.
        # If enabled but the expected feature is not actually present, reject material.
        if role in {"weak", "diagnostic", "must_have"}:
            rejected |= ~present

        # Tetracorder Ratfor excludes weak features from sumf/sumd/sumfd.
        if role != "weak":
            # This non-symettric application is faithful to
            # to the Tetracorder ratfor
            sum_fit += _zero_nan(fit * xx * weight)

            weighted_depth = depth * weight * feature_sign

            sum_depth += _zero_nan(weighted_depth)
            sum_fit_depth += _zero_nan(weighted_depth * fit)

    sum_fit[rejected] = 0.0
    sum_depth[rejected] = 0.0
    sum_fit_depth[rejected] = 0.0

    return {
        "fit": sum_fit,
        "depth": sum_depth,
        "fit_depth": sum_fit_depth,
        "rejected": rejected,
    }


def material_constraint_factor(value, values):
    if len(values) == 1:
        threshold = values[0]
        return (value >= threshold).astype(float)

    if len(values) == 2:
        reject, full = values
        return fuzzy_greater(value, (reject, full))

    raise ValueError(f"Expected 1 or 2 constraint values, got {values}")


def resolve_material_constraints(material_result, material_constraints):

    CONSTRAINT_ORDER = (
    "FITALL>",
    "DEPTHALL>",
    "FDALL>",
    "FIT>",
    "DEPTH>",
    "DEPTH-FIT>",
    "FD>",
    "FD-FIT>",
    "FD-DEPTH>",
        )


    fit = material_result["fit"].copy()
    depth = material_result["depth"].copy()
    fit_depth = material_result["fit_depth"].copy()

    constraints_by_test = {constraint["test"]: constraint for constraint in material_constraints}

    for test in CONSTRAINT_ORDER:
        # getconstraints.r:83-96 defaults every threshold to 0.0 and tp1mat
        # runs all nine blocks unconditionally, so an omitted constraint is a
        # hard ">= 0" cut. Only the three keyed on a signed quantity can bite:
        # FDALL> (tp1mat.r:751, zeroes fit, depth AND fd), FD> (873) and
        # FD-DEPTH> (901). The rest key on fit or abs(depth), never negative.
        constraint = constraints_by_test.get(test)
        values = constraint["values"] if constraint else (0.0, 0.0)
        

        if test == "FITALL>":
            factor = material_constraint_factor(fit, values)
            fit *= factor
            depth *= factor
            fit_depth *= factor

        elif test == "DEPTHALL>":
            factor = material_constraint_factor(np.abs(depth), values)
            fit *= factor
            depth *= factor
            fit_depth *= factor

        elif test == "FDALL>":
            factor = material_constraint_factor(fit_depth, values)
            fit *= factor
            depth *= factor
            fit_depth *= factor

        elif test == "FIT>":
            fit *= material_constraint_factor(fit, values)

        elif test == "DEPTH>":
            depth *= material_constraint_factor(np.abs(depth), values)

        elif test == "DEPTH-FIT>":
            depth *= material_constraint_factor(fit, values)

        elif test == "FD>":
            fit_depth *= material_constraint_factor(fit_depth, values)

        elif test == "FD-FIT>":
            fit_depth *= material_constraint_factor(fit, values)

        elif test == "FD-DEPTH>":
            fit_depth *= material_constraint_factor(depth, values)

    return {
        "fit": fit,
        "depth": depth,
        "fit_depth": fit_depth,
        "rejected": material_result["rejected"],
    }


def resolve_not(not_feature, source_result, resolved_features):
    if source_result["status"] != "valid":
        return None

    source_depth = np.abs(source_result["mod_depth"])
    source_fit = source_result["mod_fit"]

    condition = not_feature["depth_condition"]

    if condition["mode"] == "absolute":
        test_depth = source_depth

    elif condition["mode"] == "relative":
        relative_to = condition["relative_to_feature"]

        denominator = np.abs(resolved_features[relative_to]["mod_depth"])
        denominator = np.where(denominator > 1e-13, denominator, 1e-13)

        test_depth = source_depth / denominator

    else:
        raise ValueError(f"Unknown NOT depth mode: {condition['mode']}")

    return ((test_depth > condition["threshold"]) & (source_fit > not_feature["fit_threshold"]))


