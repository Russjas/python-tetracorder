"""
Class for evaluating a specific material from the tetracorder rule sets
"""
from pathlib import Path
import json
import numpy as np

from .tetracorder_ops import (GaussianConvolver,
                            linear_feature_continuum, 
                             curved_feature_continuum,
                             characterise_feature,
                             fuzzy_greater,
                             fuzzy_less,
                             continuum_test, 
                             load_references,
                             prepare_rules,
)
from .config import RRATIO_REFERENCES


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
                 mode: str = "default"
                 ):

        self.rules, self.not_definitions = self._prepare_rules(rules, mode)
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

    def evaluate(self, rule_id: str, target_spectra: np.ndarray) -> dict | None:

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
        convolvers = {}

        for rule_id, rule in self.rules.items():
            library = rule["library_records"]["SMALL"]["library"].strip("[]")
            library = LIBRARY_MAP[library]

            record = int(rule["library_records"]["SMALL"]["record"])
            reference = references[(library, record)]

            reference_wavelengths = reference["wavelengths"]

            key = (
                reference_wavelengths.shape,
                reference_wavelengths.tobytes(),
            )

            if key not in convolvers:
                convolvers[key] = GaussianConvolver(
                    lib_wl=reference_wavelengths,
                    scanner_wl=self.wavelengths,
                    scanner_fwhm=self.fwhm,
                )

            spectra[rule_id] = convolvers[key].convolve(reference["reflectance"])
        
        return spectra

    def _prepare_rratio_spectra(self, references):
        spectra = {}
        for ratio_name, key in RRATIO_REFERENCES.items():
            ratio_reference = references[key]
            convolver = GaussianConvolver(
                lib_wl=ratio_reference["wavelengths"],
                scanner_wl=self.wavelengths,
                scanner_fwhm=self.fwhm)

            spectra[ratio_name] = (convolver.convolve(ratio_reference["reflectance"]))
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

            not_definition = self.not_definitions[not_feature["source_reference"]]
            source_rule_id = not_definition["source_rule_id"]
            source_rule = self.rules[source_rule_id]
            source_rule_features = {feature["id"]: feature for feature in source_rule["features"]}
            source_feature = source_rule_features[not_feature["source_feature"]]

            source_reference_spectrum = self.rule_spectra[source_rule_id]

            source_result = resolve_feature_tests(source_reference_spectrum, source_feature,
                                                target_spectra,self.wavelengths, self.valid_bands)

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

        for rule_id in self.rules:
            results[rule_id] = self.evaluate(rule_id, target_spectra)

        return results



def fit_feature(reference, target, target_wavelengths, windows, continuum = "linear", valid_bands = None, polarity = "absorption"):
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

    Returns
    -------
    dict
        Feature fit, depth and continuum-shape measurements.

    Raises
    ------
    ValueError
        on unknown continuum type.
        
    Notes
    -----
    If a specified window is out of range of the target_wavelengths None is 
    returned as an ignore flag. This feature cannot be calculated.

    
    """
    if valid_bands is None:
        valid_bands = np.ones(target_wavelengths.shape, dtype=bool)
    else:
        valid_bands = np.asarray(valid_bands, dtype=bool)
    if valid_bands.shape != target_wavelengths.shape:
        raise ValueError("valid_bands must have the same shape as target_wavelengths")
    
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
        except ValueError as e:
            print("TARGET CONTINUUM ERROR:", e)
            return None, "invalid_data"
        
    elif continuum == "convex":
       
        for start, stop in windows:
            selected = ((target_wavelengths >= start) & (target_wavelengths <= stop))

            if not np.any(selected):
                return None, "out_of_range"

            if not np.any(selected & valid_bands):
                return None, "disabled"
        
        feature_selected = ((target_wavelengths > windows[1][1])
                            & (target_wavelengths < windows[2][0]))

        if not np.any(feature_selected):
            return None, "out_of_range"

        if not np.any(feature_selected & valid_bands):
            return None, "disabled"
        
        try:
            ref_cont = curved_feature_continuum(reference,target_wavelengths,
                                                windows,valid_bands=valid_bands)
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

    reference_cr = np.ma.asarray(ref_cont.continuum_removed)
    result["reference_area"] = float(np.ma.sum(np.ma.abs(1.0 - reference_cr)))
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

    # nvres only works over the feature and its two continuum regions.
    region = ((target_wavelengths >= left_window[0])
        & (target_wavelengths <= right_window[1]))

    wavelengths = target_wavelengths[region]
    reference = reference[region]
    ratio_reference = ratio_reference[region]
    target = np.ma.asarray(target)[..., region]
    valid = valid_bands[region]

    ratio_valid = (valid
        & np.isfinite(ratio_reference)
        & (np.abs(ratio_reference) > 1e-13))

    left = ((wavelengths >= left_window[0])
        & (wavelengths <= left_window[1])
        & ratio_valid)

    right = ((wavelengths >= right_window[0])
        & (wavelengths <= right_window[1])
        & ratio_valid)

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
        target_window = np.ma.masked_invalid(target[..., selected])
        ratio_window = ratio_reference[selected]
        ratio_data = np.broadcast_to(ratio_window, target_window.shape,)
        joint_mask = (np.ma.getmaskarray(target_window) | ~np.isfinite(ratio_data))
        target_window = np.ma.array(np.ma.getdata(target_window), mask=joint_mask,)
        ratio_window = np.ma.array(ratio_data, mask=joint_mask,)

        return (target_window.mean(axis=-1), ratio_window.mean(axis=-1),)

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

    ratioed_target = np.ma.masked_invalid(ratioed_target)

    # RRATIO-deleted channels are also unusable by bandmp.
    ratioed_target[..., ~ratio_valid] = np.ma.masked

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
    band_depth = np.ma.asarray(result["depth"]).copy()

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

    normalized_depth = np.ma.masked_invalid(normalized_depth)

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
    print("On the nvres path")
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

    fit = np.ma.asarray(result["fit"]).copy()

    depth = np.ma.asarray(result["depth"]).copy()

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
    Apply Tetracorder feature tests in Ratfor evaluation order.

    Parameters
    ----------
    feature_fit : dict
        Output of characterise_feature().

    tests : list of dict
        Parsed tests for one feature.

    Returns
    -------
    dict
        {
            "fit": ndarray,
            "depth": ndarray,
            "fit_depth": ndarray,
            "tests": dict,
        }

        fit and depth are the post-test values used by subsequent
        Tetracorder feature-role/material logic.
    """
    
    windows = feature_dict["windows"]
    continuum_type = feature_dict["continuum"]
    feat_fit, status = fit_feature(reference_spectrum, target_spectra, target_wvls, windows, continuum = continuum_type, valid_bands = valid_bands)
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

    fit = np.ma.asarray(feat_fit["fit"]).copy()
    depth = np.ma.asarray(feat_fit["depth"]).copy()

    continuum = np.ma.asarray(feat_fit["continuum"])
    left_continuum = np.ma.asarray(feat_fit["left_continuum"])
    right_continuum = np.ma.asarray(feat_fit["right_continuum"])

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

    if "r*bd>" in tests_by_name:
         # DELIBERATE DIVERGENCE from tp1mat.r (see github spectroscopy-tetracorder issue #6). 
         #The Ratfor takes the
        # reject limit from zcontlgtr(1) (getifeat.r:506/851), not zrtimesbd(1):
        #   - no lct/rct> on the feature -> reject defaults to 0.0 (811 features)
        #   - lct/rct> present           -> reject = that test's ratio limit (10 features)
        # zrtimesbd(1) is parsed (getifeat.r:1185) but never read. We honour the
        # rule's own r*bd> reject value, so this will differ from deployed Tetracorder
        # on near-threshold pixels
        value = continuum * np.abs(depth)

        factor = fuzzy_greater(value, values_for("r*bd>"))

        apply_factor("r*bd>", factor)

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

        fit = np.ma.asarray(result["mod_fit"])
        depth = np.ma.asarray(result["mod_depth"])

        feature_sign = 1.0 if result["polarity"] == "absorption" else -1.0
        weight = feature_weights[f_id]

        # Tetracorder Ratfor:
        # xx = 1.0
        # if (bdepth / xfeat <= 0.1e-5) xx = 0.0
        present = (depth / feature_sign) > 1e-6
        xx = present.astype(float)

        # featimprt > 0 means W, D or M.
        # If enabled but the expected feature is not actually present, reject material.
        if role in {"weak", "diagnostic", "must_have"}:
            rejected |= ~np.ma.filled(present, False)

        # Tetracorder Ratfor excludes weak features from sumf/sumd/sumfd.
        if role != "weak":
            # This non-symettric application is faithful to
            # to the Tetracorder ratfor
            sum_fit += np.ma.filled(fit * xx * weight, 0.0)

            weighted_depth = depth * weight * feature_sign

            sum_depth += np.ma.filled(weighted_depth, 0.0)
            sum_fit_depth += np.ma.filled(weighted_depth * fit, 0.0)

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
        if test not in constraints_by_test:
            continue

        constraint = constraints_by_test[test]
        values = constraint["values"]

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


