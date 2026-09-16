
"""
Procedural workthrough of a rule evaluation - IN PROGRESS
"""
import json
import sqlite3
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from tetracorder_ops import prepare_rules, GaussianConvolver
from time import perf_counter
from tetracorder_ops import linear_feature_continuum, curved_feature_continuum
from tetracorder_ops import characterise_feature






def load_references(db_file):
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


def get_not_reference(not_feature, references):
    LIBRARY_MAP = {"splib06": "splib06b", "sprlb06": "sprlb06b",}
    lib = not_feature["reference"]["source_library"].strip("[]")
    lib = LIBRARY_MAP[lib]

    rec = int(not_feature["reference"]["source_record"])
    key = (lib, rec)
    ref = references[key]
    return ref["reflectance"], ref["wavelengths"]


def get_reference_spectra(library_records, references):
    LIBRARY_MAP = {"splib06": "splib06b", "sprlb06": "sprlb06b",}
    lib = library_records["SMALL"]["library"].strip("[]")
    lib = LIBRARY_MAP[lib]
    rec = int(library_records["SMALL"]["record"])
    key = (lib, rec)
    ref = references[key]
    
    return ref["reflectance"], ref["wavelengths"]


def prepare_rule_references(rule, not_definitions, references, target_wavelengths, target_fwhm):
    """
    Prepare the reference spectra required by a parsed Tetracorder rule.
    
    Parameters
    ----------
    rule : dict
        Parsed, single, rule block from the JSON rules file.
    
    not_definitions : dict
        Parsed NOT-definition mapping from the JSON rules file. All NOT-definitions
    
    references : dict
        Output of load_references(); reference spectra extracted from
        the SQLite reference database.
    
    target_wavelengths : ndarray
        Wavelengths of the target spectrum or cube.
    
    target_fwhm : ndarray
        FWHM values corresponding to target_wavelengths.
    
    Returns
    -------
    main_reference : ndarray
        Main rule reference convolved onto target_wavelengths.
        Shape: (target_wavelengths.shape[0],)
    
    not_references : dict
        NOT reference spectra keyed by feature id, for example:
    
        {
            "n1a": ndarray,
            "n2a": ndarray,
        }
    
        Each spectrum is convolved onto target_wavelengths and has shape:
        (target_wavelengths.shape[0],)
    
    Notes
    -----
    1. Everything must be in microns. Reference wavelengths are stored in
       microns, so target wavelengths and FWHM must be converted before
       entering this pipeline. Reflectance must be fractional, 0-1.
    
    2. target_wavelengths.shape == target_fwhm.shape == (n_wavelengths,)
    
    3. Currently all parsed rules have a single positive reference, with
       optional NOT references. This implementation therefore assumes one
       main reference. Revisit if future rule definitions use multiple
       positive references.
    """
    main_spectrum, main_wavelengths = get_reference_spectra(
        rule["library_records"],
        references,
    )

    convolver = GaussianConvolver(
        lib_wl=main_wavelengths,
        scanner_wl=target_wavelengths,
        scanner_fwhm=target_fwhm,
    )

    main_reference = convolver.convolve(main_spectrum)

    not_references = {}

    for feature in rule["features"]:
        if feature["role"] != "not":
            continue

        not_definition = not_definitions[feature["source_reference"]]

        not_spectrum, not_wavelengths = get_not_reference(
            not_definition,
            references,)

        convolver = GaussianConvolver(
            lib_wl=not_wavelengths,
            scanner_wl=target_wavelengths,
            scanner_fwhm=target_fwhm,)
        
        convolved_spectrum = convolver.convolve(not_spectrum)
        
        source_feature = source_feature = not_definition["features"][str(feature["source_feature"])]

        not_references[feature["id"]] = {
            "spectrum": convolved_spectrum,
            "feature": source_feature,}


    return main_reference, not_references


def fit_feature(reference, target, target_wavelengths, windows, continuum = "linear"):
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
    if continuum == "linear":
        left_window = windows[0]
        right_window = windows[1]
        
        left_valid = np.any(
        (target_wavelengths >= left_window[0]) &
        (target_wavelengths <= left_window[1])
        )
        
        right_valid = np.any(
            (target_wavelengths >= right_window[0]) &
            (target_wavelengths <= right_window[1])
        )
        if not left_valid or not right_valid:
            return None, "out_of_range"
        
        
        ref_cont = linear_feature_continuum(reference, target_wavelengths,
                                            left_window,right_window,)
        target_cont = linear_feature_continuum(target, target_wavelengths,
                                               left_window, right_window,)

    elif continuum == "convex":
       
        for start, stop in windows:
            selected = ((target_wavelengths >= start)
                & (target_wavelengths <= stop))
            if not selected.any():
                return None, "out_of_range"
        try:
            ref_cont = curved_feature_continuum(reference,target_wavelengths,
                                                windows,)
        except ValueError:
            return None, "invalid_reference"
        try:
            target_cont = curved_feature_continuum(target, target_wavelengths,
                                               windows,)
        except ValueError:
            return None, "invalid_data"

    else:
        raise ValueError(f"Unknown continuum type: {continuum}")
    result = characterise_feature(ref_cont, target_cont)

    reference_cr = np.ma.asarray(ref_cont.continuum_removed)
    result["reference_area"] = float(np.ma.sum(np.ma.abs(1.0 - reference_cr)))
    return result, "valid"

def fuzzy_greater(value, thresholds):
    reject, full = thresholds
    result = np.ones_like(value, dtype=float)
    result[value <= reject] = 0.0
    transition = (value > reject) & (value < full)
    result[transition] = ((value[transition] - reject) / (full - reject))
    return result

def fuzzy_less(value, thresholds):
    reject, full = thresholds

    result = np.ones_like(value, dtype=float)

    result[value >= reject] = 0.0

    transition = (value < reject) & (value > full)
    result[transition] = ((reject - value[transition]) / (reject - full))

    return result

def continuum_test(value, limits):
    if np.isscalar(limits):
        return value >= limits

    low, high = limits
    return (value >= low) & (value <= high)

def test_ct(fit, values):
    return continuum_test(fit["continuum"], values)

def test_lct(fit, values):
    return continuum_test(fit["left_continuum"], values)

def test_rct(fit, values):
    return continuum_test(fit["right_continuum"], values)

def test_lct_rct_gt(fit, values):
    return fuzzy_greater(fit["left_right_ratio"], values)


def test_rct_lct_gt(fit, values):
    return fuzzy_greater(fit["right_left_ratio"], values)


def test_lcbbrc_gt(fit, values):
    return fuzzy_greater(fit["left_shoulder_ratio"], values)


def test_lcbbrc_lt(fit, values):
    return fuzzy_less(fit["left_shoulder_ratio"], values)


def test_rcbblc_gt(fit, values):
    return fuzzy_greater(fit["right_shoulder_ratio"], values)


def test_rcbblc_lt(fit, values):
    return fuzzy_less(fit["right_shoulder_ratio"], values)


def test_r_bd_gt(fit, values):
    return fuzzy_greater(fit["reflectance_depth"], values)

FEATURE_TESTS = {
    "ct": test_ct,
    "lct": test_lct,
    "rct": test_rct,
    
    "lct/rct>": test_lct_rct_gt,
    "rct/lct>": test_rct_lct_gt,

    "lcbbrc>": test_lcbbrc_gt,
    "lcbbrc<": test_lcbbrc_lt,
    "rcbblc>": test_rcbblc_gt,
    "rcbblc<": test_rcbblc_lt,

    "r*bd>": test_r_bd_gt,
    "weight": None,
}

def get_rule_spectrum(rule, references, target_wvls, target_fwhm):
    main_spectrum, main_wavelengths = get_reference_spectra(
        rule["library_records"],
        references,
    )

    convolver = GaussianConvolver(
        lib_wl=main_wavelengths,
        scanner_wl=target_wvls,
        scanner_fwhm=target_fwhm,
    )

    return convolver.convolve(main_spectrum)

def get_not_source_feature(not_feature, not_definitions, rules):
    not_definition = not_definitions[not_feature["source_reference"]]

    source_rule_id = not_definition["source_rule_id"]
    source_feature_number = not_feature["source_feature"]

    source_rule = next(
        rule
        for rule in rules
        if rule["id"] == source_rule_id
    )

    source_feature = next(
        feature
        for feature in source_rule["features"]
        if feature["number"] == source_feature_number
    )

    return source_rule, source_feature

def resolve_feature_tests(reference_spectrum, feature_dict, target_spectra, target_wvls):
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
    feat_fit, status = fit_feature(reference_spectrum, target_spectra, target_wvls, windows, continuum = continuum_type)
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
            if role == "must_have":
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

#%% test rules

test_rule = {'kind': 'group', 'number': 2, 'use': True, 'udata': 'reflectance', 'convolve': False, 'preratio': None, 'preprocess': None, 'algorithm': 'tricorder-primary', 'id': 'calcite.ws272.g2', 'library_records': {'SMALL': {'library': '[sprlb06]', 'record': '666'}, 'MEDIUM': {'library': '[splib06]', 'record': 'xxxx'}, 'LARGE': {'library': '[splib06]', 'record': 'xxxx'}}, 'reference_title': 'Calcite WS272                W1R1Ha', 'output_title': 'carbonate calcite WS272', 'materials': [{'name': 'calcite'}], 'identification_confidence': 8, 'features_raw': ['   f1a DLw 2.1780  2.2080  2.3770  2.4060 ct [CTHRESH4] r*bd> [RBD23]', '   n1a NOT [NOTbroadFe2]  1  0.20a  0.5   \\# NOT broad fe2+ 1-um feature', '   n2a NOT [NOTepidote]   1  0.08r1 0.3   \\# NOT epidote little 2.25 band', '\\# changed notepidote from 0.20r1 1/25/00 BWR', '\\# notepidote d/f changed from 0.15r1 0.6 1/26/00 BWR JBD', '\\#', '\\#       added broad Fe2   1/5/2000', '\\#       added not epidote 1/5/2000', '\\# Gregg Swayze determined the optional band hurts more than helps 3/95', '\\#Ow 2.1180  2.1380  2.1780  2.2080 ct 0.04', '\\# Notes:', '\\#        Feat:  1  has a weight of 0.959', '\\#        Feat:  2  has a weight of 0.041'], 'constraints': ['constraint: FD-FIT>[GLBLFDFIT] DEPTH-FIT>[GLBLDPFITg2]', 'constraint: FITALL>[GLBLFITALL]'], 'output_raw': '   output=fit depth fd\n   carbonate_calcite                     \\# Output base file name  was calcite\n   8 DN 255 = [O8DN50]\n   compress= zip', 'action': 'case 6', 'features': [{'id': 'f1a', 'number': 1, 'reference_alias': 'a', 'role': 'diagnostic', 'source_code': 'DLw', 'continuum': 'linear', 'coordinates': 'wavelength', 'windows': [[2.178, 2.208], [2.377, 2.406]], 'tests': [{'test': 'ct', 'values': 0.04}, {'test': 'r*bd>', 'values': [0.0007, 0.0011]}]}, {'id': 'n1a', 'number': 1, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTbroadFe2]', 'source_feature': 1, 'depth_condition': {'mode': 'absolute', 'threshold': 0.2}, 'fit_threshold': 0.5}, {'id': 'n2a', 'number': 2, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTepidote]', 'source_feature': 1, 'depth_condition': {'mode': 'relative', 'threshold': 0.08, 'relative_to_feature': 1}, 'fit_threshold': 0.3}]}  
o_rule = {'kind': 'group', 'number': 4, 'use': True, 'udata': 'reflectance', 'convolve': False, 'preratio': None, 'preprocess': None, 'algorithm': 'tricorder-primary', 'id': 'olivine_fo80_hs285.4b', 'library_records': {'SMALL': {'library': '[splib06]', 'record': '3696'}, 'MEDIUM': {'library': '[splib06]', 'record': 'xxxx'}, 'LARGE': {'library': '[splib06]', 'record': 'xxxx'}}, 'reference_title': 'Olivine HS285.4B Fo80 s06crj3a=b', 'output_title': 'olivine HS285.4B Fo80 =b', 'materials': [{'name': 'olivine'}], 'identification_confidence': 8, 'features_raw': ['   f1a DLw 0.680   0.777   1.720   1.850  ct .01  rct/lct> 0.6 0.7', '   n1a NOT [NOTMSNW16]    2 0.1a 0.5   \\# not melting snow 16              1.2-micron water band', '   n2a NOT [NOTH2OM5GPL]  2 0.1a 0.5   \\# not Water+Montmor SWy-2+5.01g/l  1.2-micron water band', '   n3a NOT [NOTH2OM16GPL] 2 0.1a 0.5   \\# not Water+Montmor SWy-2+16.5g/l  1.2-micron water band', '   n4a NOT [NOTWATER100]  1 0.1a 0.5   \\# not water-20C-100cm-path-model   0.74-micron water band', '\\# Notes:', '\\#         added not features to not be confused with 1.9 micron water bands t5.27b1 RNC 10/12/2022', '\\# Notes:', '\\#         Feat:  1  has a weight of', '\\#         Feat:  2  has a weight of'], 'constraints': ['constraint: FD-FIT>[GLBLFDFIT] DEPTH-FIT>[GLBLDPFITg2]', 'constraint: FITALL>[GLBLFITALL]'], 'output_raw': '   output=fit depth fd\n   olivine_fo80_hs285               \\# Output base file name\n                                    \\# previous to t5.2e1 was olivine_fo80_ki3377\n\n   8 DN 255 = 1.0000\n   compress= zip', 'action': None, 'features': [{'id': 'f1a', 'number': 1, 'reference_alias': 'a', 'role': 'diagnostic', 'source_code': 'DLw', 'continuum': 'linear', 'coordinates': 'wavelength', 'windows': [[0.68, 0.777], [1.72, 1.85]], 'tests': [{'test': 'ct', 'values': [0.01]}, {'test': 'rct/lct>', 'values': [0.6, 0.7]}]}, {'id': 'n1a', 'number': 1, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTMSNW16]', 'source_feature': 2, 'depth_condition': {'mode': 'absolute', 'threshold': 0.1}, 'fit_threshold': 0.5}, {'id': 'n2a', 'number': 2, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTH2OM5GPL]', 'source_feature': 2, 'depth_condition': {'mode': 'absolute', 'threshold': 0.1}, 'fit_threshold': 0.5}, {'id': 'n3a', 'number': 3, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTH2OM16GPL]', 'source_feature': 2, 'depth_condition': {'mode': 'absolute', 'threshold': 0.1}, 'fit_threshold': 0.5}, {'id': 'n4a', 'number': 4, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTWATER100]', 'source_feature': 1, 'depth_condition': {'mode': 'absolute', 'threshold': 0.1}, 'fit_threshold': 0.5}]}
test_rule2 ={'kind': 'group', 'number': 2, 'use': True, 'udata': 'reflectance', 'convolve': False, 'preratio': None, 'preprocess': None, 'algorithm': 'tricorder-primary', 'id': 'calcite.ws272.g2', 'library_records': {'SMALL': {'library': '[sprlb06]', 'record': '666'}, 'MEDIUM': {'library': '[splib06]', 'record': 'xxxx'}, 'LARGE': {'library': '[splib06]', 'record': 'xxxx'}}, 'reference_title': 'Calcite WS272                W1R1Ha', 'output_title': 'carbonate calcite WS272', 'materials': [{'name': 'calcite'}], 'identification_confidence': 8, 'features_raw': ['   f1a DLw 2.1780  2.2080  2.3770  2.4060 ct [CTHRESH4] r*bd> [RBD23]', '   n1a NOT [NOTbroadFe2]  1  0.20a  0.5   \\# NOT broad fe2+ 1-um feature', '   n2a NOT [NOTepidote]   1  0.08r1 0.3   \\# NOT epidote little 2.25 band', '\\# changed notepidote from 0.20r1 1/25/00 BWR', '\\# notepidote d/f changed from 0.15r1 0.6 1/26/00 BWR JBD', '\\#', '\\#       added broad Fe2   1/5/2000', '\\#       added not epidote 1/5/2000', '\\# Gregg Swayze determined the optional band hurts more than helps 3/95', '\\#Ow 2.1180  2.1380  2.1780  2.2080 ct 0.04', '\\# Notes:', '\\#        Feat:  1  has a weight of 0.959', '\\#        Feat:  2  has a weight of 0.041'], 'constraints': ['constraint: FD-FIT>[GLBLFDFIT] DEPTH-FIT>[GLBLDPFITg2]', 'constraint: FITALL>[GLBLFITALL]'], 'output_raw': '   output=fit depth fd\n   carbonate_calcite                     \\# Output base file name  was calcite\n   8 DN 255 = [O8DN50]\n   compress= zip', 'action': 'case 6', 'features': [{'id': 'f1a', 'number': 1, 'reference_alias': 'a', 'role': 'diagnostic', 'source_code': 'DLw', 'continuum': 'linear', 'coordinates': 'wavelength', 'windows': [[2.178, 2.208], [2.377, 2.406]], 'tests': [{'test': 'ct', 'values': 0.04}, {'test': 'r*bd>', 'values': [0.0007, 0.0011]}]}, {'id': 'n1a', 'number': 1, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTbroadFe2]', 'source_feature': 1, 'depth_condition': {'mode': 'absolute', 'threshold': 0.2}, 'fit_threshold': 0.5}, {'id': 'n2a', 'number': 2, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTepidote]', 'source_feature': 1, 'depth_condition': {'mode': 'relative', 'threshold': 0.08, 'relative_to_feature': 1}, 'fit_threshold': 0.3}]}
test_rule3 ={'kind': 'group', 'number': 2, 'use': True, 'udata': 'reflectance', 'convolve': False, 'preratio': None, 'preprocess': None, 'algorithm': 'tricorder-primary', 'id': 'calcite.ws272.g2', 'library_records': {'SMALL': {'library': '[sprlb06]', 'record': '666'}, 'MEDIUM': {'library': '[splib06]', 'record': 'xxxx'}, 'LARGE': {'library': '[splib06]', 'record': 'xxxx'}}, 'reference_title': 'Calcite WS272                W1R1Ha', 'output_title': 'carbonate calcite WS272', 'materials': [{'name': 'calcite'}], 'identification_confidence': 8, 'features_raw': ['   f1a DLw 2.1780  2.2080  2.3770  2.4060 ct [CTHRESH4] r*bd> [RBD23]', '   n1a NOT [NOTbroadFe2]  1  0.20a  0.5   \\# NOT broad fe2+ 1-um feature', '   n2a NOT [NOTepidote]   1  0.08r1 0.3   \\# NOT epidote little 2.25 band', '\\# changed notepidote from 0.20r1 1/25/00 BWR', '\\# notepidote d/f changed from 0.15r1 0.6 1/26/00 BWR JBD', '\\#', '\\#       added broad Fe2   1/5/2000', '\\#       added not epidote 1/5/2000', '\\# Gregg Swayze determined the optional band hurts more than helps 3/95', '\\#Ow 2.1180  2.1380  2.1780  2.2080 ct 0.04', '\\# Notes:', '\\#        Feat:  1  has a weight of 0.959', '\\#        Feat:  2  has a weight of 0.041'], 'constraints': [{'test': 'FD-FIT>', 'values': [0.3, 0.4]}, {'test': 'DEPTH-FIT>', 'values': [0.65, 0.7]}, {'test': 'FITALL>', 'values': [0.2, 0.3]}], 'output_raw': '   output=fit depth fd\n   carbonate_calcite                     \\# Output base file name  was calcite\n   8 DN 255 = [O8DN50]\n   compress= zip', 'action': 'case 6', 'features': [{'id': 'f1a', 'number': 1, 'reference_alias': 'a', 'role': 'diagnostic', 'source_code': 'DLw', 'continuum': 'linear', 'coordinates': 'wavelength', 'windows': [[2.178, 2.208], [2.377, 2.406]], 'tests': [{'test': 'ct', 'values': 0.04}, {'test': 'r*bd>', 'values': [0.0007, 0.0011]}]}, {'id': 'n1a', 'number': 1, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTbroadFe2]', 'source_feature': 1, 'depth_condition': {'mode': 'absolute', 'threshold': 0.2}, 'fit_threshold': 0.5}, {'id': 'n2a', 'number': 2, 'reference_alias': 'a', 'role': 'not', 'source_reference': '[NOTepidote]', 'source_feature': 1, 'depth_condition': {'mode': 'relative', 'threshold': 0.08, 'relative_to_feature': 1}, 'fit_threshold': 0.3}]}
#%% Test pixel setup
import spectral as sp
test_root = "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/grangegorman_84_0m00_1m00_2026-08-31_11-56-27"
test_cube = np.load(test_root + "_savgol.npy")/100.0
test_wvls = np.load(test_root + "_bands.npy")/1000.0
test_meta_path = test_root + "_metadata.json"

with open(test_meta_path, "r", encoding="utf-8") as f:
    test_meta = json.load(f)
test_fwhm = np.array([float(x) for x in test_meta["fwhm"]])[13:262]/1000.0
test_pixel = test_cube[146, 66]


#%% Startup, load files and prepare rules, nots and refs

rules_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_with_nots_parsed_normalized.json"
references_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_references.db"

with open(rules_file, "r", encoding="utf-8") as f:
    all_rules = json.load(f)


all_rules = prepare_rules(all_rules)
rules = all_rules["rules"]
not_definitions = all_rules["not_definitions"]

references = load_references(references_file)
#%%
"""
resolve positive features
        ↓
apply role semantics
        ↓
calculate material fit/depth/fd
        ↓
apply material constraints
        ↓
is material fit > 0?
        ↓ yes
evaluate NOTs
"""




#%%Rule evaluator fledgling funcs
mat_agg = {}
for rule in rules:
    #rule = test_rule #eventually this will be passed to the function this will become
    rule_features_dict = {feature["id"]: feature for feature in rule["features"]}
    rule_reference_spectrum = get_rule_spectrum(rule, references, test_wvls, test_fwhm)
    resolved_features = {}
    for f_id, feature in rule_features_dict.items():
        if feature["role"] != "not":
             result = resolve_feature_tests(rule_reference_spectrum, feature, test_cube, test_wvls)
             resolved_features[f_id] = result 
    positive_material = resolve_positive_material(rule_features_dict, resolved_features)  
    mat_agg[rule["id"]] = positive_material
#%%
from collections import Counter

summary = Counter()

for rule_id, result in mat_agg.items():
    if result is None:
        summary["no_valid_positive_features"] += 1
    elif np.all(result["rejected"]):
        summary["all_pixels_rejected"] += 1
    elif np.any(result["rejected"]):
        summary["partially_rejected"] += 1
    else:
        summary["no_pixels_rejected"] += 1

print(summary)
for rule_id, result in mat_agg.items():
    if result is None:
        continue

    for key in ("fit", "depth", "fit_depth"):
        arr = result[key]

        if not np.all(np.isfinite(arr)):
            print(rule_id, key, "contains non-finite values")
            
for rule in rules:
    rule_features_dict = {feature["id"]: feature for feature in rule["features"]}
    resolved = {}

    rule_reference_spectrum = get_rule_spectrum(rule, references, test_wvls, test_fwhm)

    for f_id, feature in rule_features_dict.items():
        if feature["role"] != "not":
            resolved[f_id] = resolve_feature_tests(rule_reference_spectrum, feature, test_cube, test_wvls)

    weights = resolve_feature_weights(rule_features_dict, resolved)

    total = sum(weights.values())

    if weights and not np.isclose(total, 1.0) and total != 0.0:
        print(rule["id"], total, weights)
#%%
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

mat_agg = {}

for rule in rules:
    rule_features_dict = {feature["id"]: feature for feature in rule["features"]}
    rule_reference_spectrum = get_rule_spectrum(rule, references, test_wvls, test_fwhm)

    resolved_features = {}

    for f_id, feature in rule_features_dict.items():
        if feature["role"] != "not":
            resolved_features[f_id] = resolve_feature_tests(rule_reference_spectrum, feature, test_cube, test_wvls)

    positive_material = resolve_positive_material(rule_features_dict, resolved_features)
    mat_agg[rule["id"]] = positive_material

profiler.disable()

stats = pstats.Stats(profiler)
stats.sort_stats("cumulative")
stats.print_stats(30)

#%%
for key, im in positive_material.items():
    plt.figure()
    plt.imshow(im)
    plt.title(key)     
#%%

rej = positive_material["rejected"]

print("Rejected pixels:", np.sum(rej))
print("Total pixels:", rej.size)
print("Rejected fraction:", np.mean(rej))

print("fit nonzero in rejected:", np.count_nonzero(positive_material["fit"][rej]))
print("depth nonzero in rejected:", np.count_nonzero(positive_material["depth"][rej]))
print("fit_depth nonzero in rejected:", np.count_nonzero(positive_material["fit_depth"][rej])) 
feature_weights = resolve_feature_weights(rule_features_dict, resolved_features)

for f_id, weight in feature_weights.items():
    feature = rule_features_dict[f_id]
    result = resolved_features[f_id]

    print(
        f_id,
        feature["role"],
        result["status"],
        "area =", result.get("reference_area"),
        "weight =", weight,
    )

print("sum =", sum(feature_weights.values()))       

for f_id, feature in rule_features_dict.items():
    if feature["role"] == "not":
        continue

    result = resolved_features[f_id]

    if result["status"] != "valid":
        print(f_id, feature["role"], result["status"])
        continue

    depth = np.ma.asarray(result["mod_depth"])
    feature_sign = 1.0 if result["polarity"] == "absorption" else -1.0
    present = (depth / feature_sign) > 1e-6

    print(f_id, feature["role"], "present fraction =", np.mean(np.ma.filled(present, False)))
#%% Not work, to revist when I get there is resolution order
for f_id, feature in rule_features_dict.items():
    if feature["role"] == "not":
        source_rule, source_feature = get_not_source_feature(feature, not_definitions, rules)
        not_reference_spectrum = get_rule_spectrum(source_rule, references, test_wvls, test_fwhm)
        source_result = resolve_feature_tests(not_reference_spectrum, source_feature, test_cube, test_wvls)
        result = resolve_not(feature, source_result, resolved_features)
        resolved_features[f_id] = result 

#%%
for f_id, feature in rule_features_dict.items():
    if feature["role"] == "not":
        print(f_id, feature["depth_condition"])


#%%
for rule in rules:
    #rule = test_rule3 #eventually this will be passed to the function this will become
    rule_features_dict = {feature["id"]: feature for feature in rule["features"]}
    rule_reference_spectrum = get_rule_spectrum(rule, references, test_wvls, test_fwhm)
    resolved_features = {}
    for f_id, feature in rule_features_dict.items():
        if feature["role"] != "not":
             result = resolve_feature_tests(rule_reference_spectrum, feature, test_cube, test_wvls)
        else:
            source_rule, source_feature = get_not_source_feature(feature, not_definitions, rules)
            not_reference_spectrum = get_rule_spectrum(source_rule, references, test_wvls, test_fwhm)
            result = resolve_feature_tests(not_reference_spectrum, source_feature, test_cube, test_wvls)
        resolved_features[f_id] = result    
        
#%%
for rule in rules:
    features_by_number = {feature["number"]: feature["id"] for feature in rule["features"]}

    for feature in rule["features"]:
        if feature["role"] != "not":
            continue

        condition = feature["depth_condition"]

        if condition["mode"] == "relative":
            feature_number = condition["relative_to_feature"]
            condition["relative_to_feature"] = features_by_number[feature_number]

