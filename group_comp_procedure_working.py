"""
Messing space from group evaluation dev


"""
import json

import spectral as sp
import numpy as np

from material_evaluation import MaterialEvaluator

#%% test data set up


test_root = "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/grangegorman_84_0m00_1m00_2026-08-31_11-56-27"
test_cube = np.load(test_root + "_savgol.npy")/100.0
test_wvls = np.load(test_root + "_bands.npy")/1000.0
test_meta_path = test_root + "_metadata.json"

with open(test_meta_path, "r", encoding="utf-8") as f:
    test_meta = json.load(f)
test_fwhm = np.array([float(x) for x in test_meta["fwhm"]])[13:262]/1000.0
test_pixel = test_cube[146, 66]


#%% test rule ids
RULE_IDS = ['calcite.ws272.g2',
            'olivine_fo80_hs285.4b',
            "clinochlore.fe.gds157",
            "pyroxene.hypersthene.pyx02.e.g1" 
    ]

#%% Startup, load files init evaluator
from time import perf_counter

t0 = perf_counter()
rules_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_as_dict.json"
references_file = "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/tetracorder_rules_references.db"
evaluator = MaterialEvaluator(
    rules_file,
    references_file,
    test_wvls,
    test_fwhm,
    mode="default")


print(f"Evaluator initialisation: {perf_counter() - t0:.3f} s")
#%% map
"""
ALL ORDINARY MATERIALS
        ↓
your current per-material evaluator
        ↓
fit/depth/fd after NOTs
        ↓
GROUP COMPETITION 
        ↓
highest fit wins each pixel in each group
        ↓
optional class override # Not implemented, no rules use
        ↓
zero all losing materials #? probably unnecessary in this implementation
        ↓
winning material has action?<-current
        ↓ yes
run CASE(s)
        ↓
case-specific material evaluation / resolution
        ↓
final outputs
"""
#%% Evaluate material and assemble into groups
materials = evaluator.rules
group0 = {}
groups = {}

for mid, rule in evaluator.rules.items():
    if rule["kind"] != "group":
        continue

    result = evaluator.evaluate(mid, test_cube)

    if rule["number"] == 0:
        group0[mid] = result
    else:
        groups.setdefault(rule["number"], {})[mid] = result
    
#%%
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

group_winners = {}

for group_num, group_results in groups.items():
    group_winners[group_num] = resolve_group(
        group_results,
        group0=group0,
    )
    
#%%

for group_num, result in group_winners.items():

    if result is None:
        continue

    winner = result["winner"]

    for mid in np.unique(winner[winner != None]):
        action = evaluator.rules[mid].get("action")

        if action is None:
            continue

        else:
            print(action)
        winner_mask = winner == mid

        # action e.g. "case 6"
        # run that case only for winner_mask









