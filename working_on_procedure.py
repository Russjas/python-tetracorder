
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
        "linear" or "curved". Continuum calculation type specified in the rule

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
        return None
    
    if continuum == "linear":
        ref_cont = linear_feature_continuum(reference, target_wavelengths,
                                            left_window,right_window,)
        target_cont = linear_feature_continuum(target, target_wavelengths,
                                               left_window, right_window,)

    elif continuum == "curved":
        ref_cont = curved_feature_continuum(reference,target_wavelengths,
                                            left_window, right_window,)
        target_cont = curved_feature_continuum(target, target_wavelengths,
                                               left_window, right_window,)

    else:
        raise ValueError(f"Unknown continuum type: {continuum}")

    return characterise_feature(ref_cont, target_cont,)

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

def test_ct(feature_fit, threshold):
    return feature_fit["continuum"] >= threshold


def test_lct(feature_fit, threshold):
    return feature_fit["left_continuum"] >= threshold


def test_rct(feature_fit, threshold):
    return feature_fit["right_continuum"] >= threshold

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

rules_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_with_nots_parsed.json"
references_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_references.db"

with open(rules_file, "r", encoding="utf-8") as f:
    all_rules = json.load(f)

# =============================================================================
# import re
# 
# CONSTRAINT_RE = re.compile(
#     r"(FITALL>|FIT>|DEPTHALL>|DEPTH-FIT>|FD-FIT>|FDALL>|FD-DEPTH>|DEPTH>|FD>)"
#     r"\s*"
#     r"(\[[^\]]+\]|[-+]?\d*\.?\d+(?:\s+[-+]?\d*\.?\d+)?)"
# )
# 
# def parse_constraints(constraints):
#     parsed = []
# 
#     for line in constraints:
#         line = line.removeprefix("constraint:").strip()
# 
#         for test, values in CONSTRAINT_RE.findall(line):
#             parsed.append({
#                 "test": test,
#                 "values": values.strip(),
#             })
# 
#     return parsed
# 
# for rule in all_rules["rules"]:
#         if rule.get("constraints"):
#             rule["constraints"] = parse_constraints(rule["constraints"])
# 
# output_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_with_nots_parsed.json"
# 
# =============================================================================
# =============================================================================
# with open(output_file, "w", encoding="utf-8") as f:
#     json.dump(all_rules, f, indent=2)
# =============================================================================
    
#%%

all_rules = prepare_rules(all_rules)
rules = all_rules["rules"]
nots = all_rules["not_definitions"]

references = load_references(references_file)

for rule in rules:
    if rule["id"] == 'calcite.ws272.g2':
        print(rule)


#%%Rule evaluator fledgling funcs

rule = test_rule3 #eventually this will be passed to the function this will become
rule_features_dict = {feature["id"]: feature for feature in rule["features"]}
main_reference, not_references = prepare_rule_references(rule, nots, references, test_wvls, test_fwhm)

fits = {}
for feature_id, feature in rule_features_dict.items():
    #(reference, target, target_wavelengths, windows, continuum = "linear")
    
    if feature["role"] == "not":
        reference = not_references[feature_id]["spectrum"] 
        windows = not_references[feature_id]["feature"]["windows"]
        continuum = not_references[feature_id]["feature"]["continuum"]
        feat_fit = fit_feature(reference, test_cube, test_wvls, windows, continuum = continuum)
    else:
        windows = feature["windows"]
        continuum = feature["continuum"]
        feat_fit = fit_feature(main_reference, test_cube, test_wvls, windows, continuum = continuum)
    fits[feature_id] = feat_fit

test_results = {}

for feature_id in fits.keys():
    if fits[feature_id] is None:
        continue

    tests = rule_features_dict[feature_id].get("tests")
    if tests is None:
        continue

    test_results[feature_id] = {}

    for test in tests:
        check = test["test"]
        values = test["values"]

        if isinstance(values, (list, tuple)) and len(values) == 1:
            values = values[0]

        test_func = FEATURE_TESTS[check]
        if test_func is not None:
            result = test_func(fits[feature_id], values)
        else:
            result = None

        test_results[feature_id][check] = result

#%%

