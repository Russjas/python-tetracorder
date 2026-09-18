"""
Top level module to evaluate sample spectra against the full Tetracorder numerical pipeline
"""
from pathlib import Path
import numpy as np

from .material_evaluation import MaterialEvaluator 
from .config import NOGROUP0 




class GroupEvaluator:
    def __init__(self, target_spectra, target_wavelengths, target_fwhm, mode = "default", target_valid_bands = None, reference_file = None, rules_file = None ):
        self.target =  target_spectra
        self.wavelengths = target_wavelengths
        self.fwhm = target_fwhm
        self.mode = mode
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
                mode=mode)
        self.group_winners, self.case_winners = self.evaluate()

    def evaluate(self):

        # Evaluate all of the materials that have rules provided. 
        # Any material exclusion rule logic could go here        
        group0 = {}
        groups = {}

        for mid, rule in self.evaluator.rules.items():
            if rule["kind"] != "group":
                continue
            print(f"Evaluating {mid}")
            result = self.evaluator.evaluate(mid, self.target)
            if rule["number"] == 0:
                group0[mid] = result
            else:
                groups.setdefault(rule["number"], {})[mid] = result

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
            #search for case definitions in the rule file
            if rule["kind"] != "case":
                continue
            case_num = int(rule["number"])
            if case_num not in case_masks:
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