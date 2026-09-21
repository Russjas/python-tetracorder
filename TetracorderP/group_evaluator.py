"""
Top level module to evaluate sample spectra against the full Tetracorder numerical pipeline
"""
from pathlib import Path
import gzip
import re
import json

import numpy as np

from .material_evaluation import MaterialEvaluator 
from .config import NOGROUP0, VARIABLE_PRESETS, OUTPUT_DIRS




class GroupEvaluator:
    def __init__(self, target_spectra, target_wavelengths, target_fwhm, mode = "default", target_valid_bands = None, reference_file = None, rules_file = None,
                 disabled_groups = None, disabled_cases = None, disabled_materials = None,
                 temperature = None, pressure = None ):
        self.temperature = temperature    # (min, max) Kelvin of the data set, or None
        self.pressure = pressure          # (min, max) bar of the data set, or None
        self.target =  target_spectra
        self.wavelengths = target_wavelengths
        self.fwhm = target_fwhm
        self.mode = mode

        self.disabled_groups = {int(g) for g in (disabled_groups or ())}
        self.disabled_cases = {int(c) for c in (disabled_cases or ())}
        self.disabled_materials = set(disabled_materials or ())

        self.valid_bands = target_valid_bands
        if self.valid_bands is None:
            self.valid_bands = np.ones(self.wavelengths.shape, dtype=bool)
        self.reference_file = reference_file
        self.rules_file = rules_file
        self._find_files()
        self.evaluator = MaterialEvaluator(
                self.rules_file,
                self.reference_file,
                self.wavelengths,
                self.fwhm,
                target_valid_bands = self.valid_bands,
                mode=mode,
                disabled_materials = self.disabled_materials)
        self.disabled_materials |= self._resolve_physical_and_disable()
        self.evaluator.disabled_materials = self.disabled_materials
        self.group_winners, self.case_winners = self.evaluate()

    def _resolve_physical_and_disable(self):
        """Materials whose declared temperature / pressure range excludes the data.

        applygtpconstraints.r disables a material outright when the data range
        lies entirely outside [limit1, limit4]; the inner two limits are parsed
        by the Ratfor and never read. Limits and data ranges are Kelvin and bar.
        A material with no declared limit, or a run with no declared condition,
        is left enabled.
        """
        disabled = set()
        conditions = (("temperature", self.temperature), ("pressure", self.pressure))

        for mid, rule in self.evaluator.rules.items():
            if mid in self.disabled_materials:
                continue
            for name, data_range in conditions:
                limits = rule.get("physical", {}).get(name)
                if not limits or data_range is None:
                    continue
                low, high = limits[0], limits[3]
                if low is not None and data_range[1] < low:
                    disabled.add(mid)
                    break
                if high is not None and data_range[0] > high:
                    disabled.add(mid)
                    break
        return disabled

    def evaluate(self):

        # Evaluate all of the materials that have rules provided. 
        # Any material exclusion rule logic could go here        
        group0 = {}
        groups = {}
        count = 1
        for mid, rule in self.evaluator.rules.items():
            if rule["kind"] != "group":
                continue
            if int(rule["number"]) in self.disabled_groups:
                continue
            if mid in self.disabled_materials:
                continue
            print(f"Evaluating {count} of {len(self.evaluator.rules)}: {mid}")
            result = self.evaluator.evaluate(mid, self.target)
            if rule["number"] == 0:
                group0[mid] = result
            else:
                groups.setdefault(rule["number"], {})[mid] = result
            count +=1

        # Evaluate each group to assemble the group winners
        group_winners = {}
        
        for group_num, group_results in groups.items():
            include_group0 = group_num not in NOGROUP0
            group_winners[group_num] = resolve_group(group_results, group0=group0 if include_group0 else None,)

        # Evaluate any cases held by the group winners        
        case_masks = {}
        
        for group_num, result in group_winners.items():
            # Ignore groups that did not have a valid winner
            if result is None:
                continue
        
            winner = result["winner"]
        
            for mid in np.unique(winner[winner != None]):
                #For each unqique winner, check if it has an action
                action = self.evaluator.rules[mid].get("action")
                if action is None:
                    continue
                #if there is an action, parse it
                parts = action.split()
                if parts[0] != "case":
                    #TODO: implement non case actions (except play sound, because no)
                    continue

                # In the ratfor a case is only evaluated when: group winner is valid AND has a non-zero depth
                winner_mask = (winner == mid) & (np.abs(result["depth"]) > 0.0)

                # most case actions have the form action: case 1
                # some have the for action: case 1 2 3 4 5
                # so we iterate to resolve every case
                for case_num in [int(p) for p in parts[1:]]:
                    if case_num not in case_masks:
                        case_masks[case_num] = winner_mask.copy()
                    else:
                        case_masks[case_num] |= winner_mask

        # Evaluate relevant cases
        cases = {}
        
        for mid, rule in self.evaluator.rules.items():
            if rule["kind"] != "case":
                continue
            case_num = int(rule["number"])
            if case_num not in case_masks:
                continue
            if case_num in self.disabled_cases or mid in self.disabled_materials:
                continue
            # evaluate those cases as if they were group rules
            # and store the results
            cases.setdefault(case_num, {})[mid] = self.evaluator.evaluate(mid, self.target)

        #evaluate cases to find case winners
        case_winners = {}
        
        for case_num, case_results in cases.items():
            case_winners[case_num] = resolve_case(case_results, case_masks[case_num])
        return group_winners, case_winners


    def _find_files(self):
        """
        Find missing Tetracorder support files.

        Explicitly supplied file paths are preserved. Only missing files
        are searched for.
        """
        if self.reference_file is None:
            self.reference_file = self._find_file(
                "tetracorder_rules_references.db")

        if self.rules_file is None:
            self.rules_file = self._find_file("tetracorder_rules_as_dict.json")


    def _find_file(self, filename):
        """
        Search from the current module root.
        """
        module_root = Path(__file__).resolve().parent.parent

        for path in module_root.rglob(filename):
            if path.is_file():
                return str(path)

        raise FileNotFoundError(f"Could not find required Tetracorder file: {filename}")

    #==========writing functions =======================================================
    def write_like_tetracorder(self, output_dir):
        """
        Write every enabled material's fit, depth and fd images as Tetracorder
        does: <dir>/<name>.<plane>.gz (VICAR label + 8-bit image, gzipped) and
        <name>.<plane>.gz.hdr. Pixels a material did not win are 0. Group-0
        materials are written into every group directory that includes group 0.
        """
        output_dir = Path(output_dir)
        rules = self.evaluator.rules
        nl, ns = self.target.shape[:2]
        lblsiz = ns if ns >= 299 else ns * (299 // ns + 1)      # creatoutfiles.r
        label = (f"LBLSIZE={lblsiz}  FORMAT='BYTE'  TYPE='IMAGE'  RECSIZE={ns}  "
                f"ORG='BSQ'  NL={nl}  NS={ns}  NB=1  ").encode().ljust(lblsiz)

        # groups that write output, and those that take group 0 (cubecorder.r:448-468)
        live_groups = {int(r["number"]) for r in rules.values() if r["kind"] == "group"} \
            - {0} - set(self.disabled_groups)
        takes_group0 = sorted(live_groups - NOGROUP0)

        for rid, rule in rules.items():
            if rid in self.disabled_materials:
                continue
            number = int(rule["number"])
            if rule["kind"] == "case":
                if number in self.disabled_cases:
                    continue
                targets = [(OUTPUT_DIRS["case"][number], self.case_winners.get(number))]
            elif number == 0:
                targets = [(OUTPUT_DIRS["group"][g], self.group_winners.get(g))
                        for g in takes_group0]
            elif number in live_groups:
                targets = [(OUTPUT_DIRS["group"][number], self.group_winners.get(number))]
            else:
                continue

            # output block: "output=fit depth fd" / "<name>" / "8 DN 255 = <value>"
            lines = [l.split("\\#")[0].strip() for l in rule["output_raw"].splitlines()]
            lines = [l for l in lines if l]
            name = lines[1].split()[0]
            dn, value = re.match(r"\d+\s+DN\s+(\d+)\s*=\s*(\S+)", lines[2]).groups()
            value = VARIABLE_PRESETS[self.mode].get(value, value)
            depth_scale = int(dn) / float(value)

            for subdir, result in targets:
                (output_dir / subdir).mkdir(parents=True, exist_ok=True)
                won = (np.asarray(result["winner"], dtype=object) == rid
                    if result is not None else np.zeros((nl, ns), dtype=bool))

                for plane, key, scale, suffix in (("fit", "fit", 255.0, "FIT"),
                                                ("depth", "depth", depth_scale, "DEPTHS"),
                                                ("fd", "fit_depth", depth_scale, "FIT*DEPTH")):
                    values = (np.ma.filled(np.ma.asarray(result[key], dtype=np.float32), 0.0)
                            if result is not None else np.zeros((nl, ns), np.float32))
                    x = np.where(won, values, 0.0).astype(np.float32) * np.float32(scale)
                    image = np.clip(np.floor(np.nan_to_num(x).astype(np.float64) + 0.5), 0, 255)

                    path = output_dir / subdir / f"{name}.{plane}"
                    with gzip.GzipFile(f"{path}.gz", "wb", compresslevel=6, mtime=0) as f:
                        f.write(label + image.astype(np.uint8).tobytes())
                    Path(f"{path}.gz.hdr").write_text(
                        f"ENVI\ndescription = {{\n  {rule['output_title']} {suffix}\n  }}\n"
                        f"samples = {ns}\nlines   = {nl}\nbands   = 1\n"
                        f"header offset = {lblsiz}\nfile compression = 1\n"
                        "file type = ENVI Standard\ndata type = 1\ninterleave = bsq\n"
                        "sensor type = spectral data\nbyte order = 0\n"
                        "wavelength units = Micrometers\n")

    def write_npz(self, path):
        """
        Save every group and case result to one compressed .npz:

            materials            str, index -> material id; 0 = "" (nothing detected)
            group_N_material     int16 code into materials
            group_N_fit          float32
            group_N_depth        float32
            group_N_fit_depth    float32
            case_N_...           the same for each case
            wavelengths, valid_bands
            meta                 JSON string: mode, files, conditions, disabled lists

        Read back with np.load(path); no pickling is needed or allowed.
        """
        materials = [""] + list(self.evaluator.rules)
        code = {mid: i for i, mid in enumerate(materials)}
        arrays = {"materials": np.array(materials, dtype=str),
                "wavelengths": np.asarray(self.wavelengths, dtype=np.float32),
                "valid_bands": np.asarray(self.valid_bands, dtype=bool)}

        for kind, results in (("group", self.group_winners), ("case", self.case_winners)):
            for number, result in results.items():
                if result is None:
                    continue
                winner = np.asarray(result["winner"], dtype=object)
                codes = np.zeros(winner.shape, dtype=np.int16)
                for mid in set(winner[winner != None]):  # noqa: E711
                    codes[winner == mid] = code[mid]
                arrays[f"{kind}_{number}_material"] = codes
                for key in ("fit", "depth", "fit_depth"):
                    arrays[f"{kind}_{number}_{key}"] = np.asarray(
                        np.ma.filled(np.ma.asarray(result[key]), 0.0), dtype=np.float32)

        arrays["meta"] = np.array(json.dumps({
            "mode": self.mode,
            "rules_file": str(self.rules_file),
            "reference_file": str(self.reference_file),
            "temperature_K": list(self.temperature) if self.temperature is not None else None,
            "pressure_bar": list(self.pressure) if self.pressure is not None else None,
            "disabled_groups": sorted(self.disabled_groups),
            "disabled_cases": sorted(self.disabled_cases),
            "disabled_materials": sorted(self.disabled_materials),
        }))
        np.savez_compressed(Path(path), **arrays)


    


def resolve_group(group_results, group0=None):

    candidates = {mid: result for mid, result in group_results.items()
        if result is not None}

    if group0 is not None:
        candidates.update({mid: result for mid, result in group0.items()
            if result is not None})

    if not candidates:
        return None
    
    material_ids = list(candidates)

    fit_stack = np.stack([candidates[mid]["fit"] for mid in material_ids])
    depth_stack = np.stack([candidates[mid]["depth"] for mid in material_ids])
    fd_stack = np.stack([candidates[mid]["fit_depth"] for mid in material_ids])


    #find winning material for the group
    winner_index = np.argmax(fit_stack, axis=0)
    # Tetracorder has a complicated 2nd best winner with a class overide system
    # implemented here. But none of the existing rules has a constraint: class
    # so it is never used. Not implemented here until needed   

    # find the CORRESPONDING fit, depth and depth*fit to the winner
    winner_fit = np.take_along_axis(fit_stack, winner_index[None, ...], axis=0)[0]
    winner_depth = np.take_along_axis(depth_stack, winner_index[None, ...], axis=0)[0]
    winner_fd = np.take_along_axis(fd_stack, winner_index[None, ...], axis=0)[0]

    # Array of winner names
    material_ids = np.asarray(material_ids, dtype=object)
    winner = material_ids[winner_index]

    # Mask of pixels where no material in this group was valid
    # this avoids assigning the pixel to the first listed material in the group
    detected = winner_fit > 0.0

    winner = np.where(detected, winner, None)
    winner_fit = np.where(detected, winner_fit, 0.0)
    winner_depth = np.where(detected, winner_depth, 0.0)
    winner_fd = np.where(detected, winner_fd, 0.0)

    return {
        "winner": winner,
        "fit": winner_fit,
        "depth": winner_depth,
        "fit_depth": winner_fd,
    }

def resolve_case(case_results, active_mask):

    candidates = {mid: result for mid, result in case_results.items()
        if result is not None}

    if not candidates:
        return None

    material_ids = list(candidates)

    fit_stack = np.stack([candidates[mid]["fit"] for mid in material_ids])
    depth_stack = np.stack([candidates[mid]["depth"] for mid in material_ids])
    fd_stack = np.stack([candidates[mid]["fit_depth"] for mid in material_ids])

    # Highest-fit case material wins each pixel
    winner_index = np.argmax(fit_stack, axis=0)

    winner_fit = np.take_along_axis(fit_stack, winner_index[None, ...], axis=0)[0]
    winner_depth = np.take_along_axis(depth_stack, winner_index[None, ...], axis=0)[0]
    winner_fd = np.take_along_axis(fd_stack, winner_index[None, ...], axis=0)[0]

    material_ids = np.asarray(material_ids, dtype=object)
    winner = material_ids[winner_index]

    # A case result is valid only where:
    # 1. the case was invoked by a group winner
    # 2. at least one case material survived evaluation
    detected = active_mask & (winner_fit > 0.0)

    winner = np.where(detected, winner, None)
    winner_fit = np.where(detected, winner_fit, 0.0)
    winner_depth = np.where(detected, winner_depth, 0.0)
    winner_fd = np.where(detected, winner_fd, 0.0)

    return {
        "winner": winner,
        "fit": winner_fit,
        "depth": winner_depth,
        "fit_depth": winner_fd,
    }