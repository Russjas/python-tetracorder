"""
Messing space from group evaluation dev


"""
import json

import spectral as sp
import numpy as np


from src.group_evaluator import GroupEvaluator


root_name = "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/dh-64-03-carroweagh_0_0m00_1m00_2026-08-31_11-50-53"


test_cube = np.load(root_name + "_savgol.npy")/100.0
test_wvls = np.load(root_name + "_bands.npy")/1000.0
test_meta_path = root_name + "_metadata.json"

with open(test_meta_path, "r", encoding="utf-8") as f:
    test_meta = json.load(f)
test_fwhm = np.array([float(x) for x in test_meta["fwhm"]])[13:262]/1000.0
test_pixel = test_cube[146, 66]


tetra = GroupEvaluator(test_cube, test_wvls, test_fwhm)
#%%
rules = tetra.evaluator.rules
for rid, rule in rules.items():
    if rule["algorithm"]!= "tricorder-primary":
        print(rid)
#%%
print(rules["red.edge.shift.2"])
#%% test data set up
tests = {"test1": "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/18-kn-03_31_0m00_1m00_2026-08-31_11-32-06",
"test2": "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/baltinglas-granint_2_1m00_2m00_2026-08-31_11-39-08",
"test3": "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/dh-64-03-carroweagh_0_0m00_1m00_2026-08-31_11-50-53",
"test4": "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/dh-302-silvermines_2_0m00_1m00_2026-08-31_11-45-38",}
for test_name, test_root in tests.items():
#test_name = "test1"
#test_root = "C:/Users/Hyperspectral/Documents/HS_Data/Exhibit_boxes/20260831/SWIR/Processed_exhibits/grangegorman_84_0m00_1m00_2026-08-31_11-56-27"
    test_cube = np.load(test_root + "_savgol.npy")/100.0
    test_wvls = np.load(test_root + "_bands.npy")/1000.0
    test_meta_path = test_root + "_metadata.json"
    
    with open(test_meta_path, "r", encoding="utf-8") as f:
        test_meta = json.load(f)
    test_fwhm = np.array([float(x) for x in test_meta["fwhm"]])[13:262]/1000.0
    test_pixel = test_cube[146, 66]
    
    
    #% Startup, load files init evaluator
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
    
    #% Evaluate material and assemble into groups
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
    group_winners = {}
    
    for group_num, group_results in groups.items():
        group_winners[group_num] = resolve_group(
            group_results,
            group0=group0,
        )
    
    
    case_masks = {}
    
    for group_num, result in group_winners.items():
    
        if result is None:
            continue
    
        winner = result["winner"]
    
        for mid in np.unique(winner[winner != None]):
    
            action = evaluator.rules[mid].get("action")
    
            if action is None:
                continue
    
            _, case_num = action.split()
            case_num = int(case_num)
    
            winner_mask = (winner == mid) & (np.abs(result["depth"]) > 0.0)
    
            if case_num not in case_masks:
                case_masks[case_num] = winner_mask.copy()
            else:
                case_masks[case_num] |= winner_mask
    cases = {}
    
    for mid, rule in evaluator.rules.items():
    
        if rule["kind"] != "case":
            continue
    
        case_num = rule["number"]
    
        if case_num not in case_masks:
            continue
    
        cases.setdefault(case_num, {})[mid] = evaluator.evaluate(mid, test_cube)
    
    case_winners = {}
    
    for case_num, case_results in cases.items():
        case_winners[case_num] = resolve_case(
            case_results,
            case_masks[case_num],
        )
    
    
    def build_final_material_outputs(group_winners, case_winners):
    
        final_outputs = {}
    
        for group_num, result in group_winners.items():
    
            if result is None:
                continue
    
            winner = result["winner"]
    
            for mid in np.unique(winner[winner != None]):
    
                mask = winner == mid
    
                if mid not in final_outputs:
                    final_outputs[mid] = {
                        "fit": np.zeros_like(result["fit"]),
                        "depth": np.zeros_like(result["depth"]),
                        "fit_depth": np.zeros_like(result["fit_depth"]),
                    }
    
                final_outputs[mid]["fit"][mask] = result["fit"][mask]
                final_outputs[mid]["depth"][mask] = result["depth"][mask]
                final_outputs[mid]["fit_depth"][mask] = result["fit_depth"][mask]
    
        for case_num, result in case_winners.items():
    
            if result is None:
                continue
    
            winner = result["winner"]
    
            for mid in np.unique(winner[winner != None]):
    
                mask = winner == mid
    
                if mid not in final_outputs:
                    final_outputs[mid] = {
                        "fit": np.zeros_like(result["fit"]),
                        "depth": np.zeros_like(result["depth"]),
                        "fit_depth": np.zeros_like(result["fit_depth"]),
                    }
    
                final_outputs[mid]["fit"][mask] = result["fit"][mask]
                final_outputs[mid]["depth"][mask] = result["depth"][mask]
                final_outputs[mid]["fit_depth"][mask] = result["fit_depth"][mask]
    
        return final_outputs
    results = build_final_material_outputs(
        group_winners,
        case_winners,
    )
    #%
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Patch
    
    
    group_num = 2
    fit_threshold = 0.90
    top_n = 20
    
    
    # Resolve both versions
    g2_without_group0 = resolve_group(groups[group_num])
    g2_with_group0 = resolve_group(groups[group_num], group0=group0)
    
    winner_no_g0 = g2_without_group0["winner"]
    winner_with_g0 = g2_with_group0["winner"]
    
    
    # Common material list so colours mean the same thing in both winner maps
    materials = np.unique(np.concatenate([
        winner_no_g0[winner_no_g0 != None],
        winner_with_g0[winner_with_g0 != None]
    ]))
    
    material_index = {mid: i for i, mid in enumerate(materials)}
    
    
    def make_index_map(winner):
    
        index_map = np.full(winner.shape, -1, dtype=int)
    
        for mid, i in material_index.items():
            index_map[winner == mid] = i
    
        return np.ma.masked_where(index_map < 0, index_map)
    
    
    map_no_g0 = make_index_map(winner_no_g0)
    map_with_g0 = make_index_map(winner_with_g0)
    
    
    # Winner frequencies, used to select legend entries
    vals, counts = np.unique(
        winner_no_g0[winner_no_g0 != None],
        return_counts=True
    )
    
    order = np.argsort(counts)[::-1]
    top_materials = vals[order][:top_n]
    top_counts = counts[order][:top_n]
    
    
    # Categorical colour map
    cmap = plt.get_cmap("tab20", len(materials))
    
    
    # Threshold the winning fit
    fit = g2_without_group0["fit"]
    thresholded_fit = np.ma.masked_where(fit < fit_threshold, fit)
    
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 10))
    
    
    axes[0].imshow(
        map_no_g0,
        cmap=cmap,
        vmin=0,
        vmax=len(materials) - 1,
        interpolation="nearest"
    )
    axes[0].set_title(f"Group {group_num} without group 0")
    axes[0].axis("off")
    
    
    axes[1].imshow(
        map_with_g0,
        cmap=cmap,
        vmin=0,
        vmax=len(materials) - 1,
        interpolation="nearest"
    )
    axes[1].set_title(f"Group {group_num} with group 0")
    axes[1].axis("off")
    
    
    im = axes[2].imshow(
        thresholded_fit,
        vmin=fit_threshold,
        vmax=1.0,
        interpolation="nearest"
    )
    axes[2].set_title(f"Winning fit ≥ {fit_threshold}")
    axes[2].axis("off")
    
    
    cbar = fig.colorbar(
        im,
        ax=axes[2],
        fraction=0.046,
        pad=0.04
    )
    cbar.set_label("Fit")
    
    
    # Top-N winner legend
    handles = [
        Patch(
            facecolor=cmap(material_index[mid]),
            label=f"{mid} ({count})"
        )
        for mid, count in zip(top_materials, top_counts)
    ]
    
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        fontsize=8
    )
    
    
    plt.tight_layout(rect=[0, 0, 0.87, 1])
    plt.savefig(f"{test_name}.jpg")
    plt.show()
