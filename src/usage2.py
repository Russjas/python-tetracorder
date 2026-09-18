import numpy as np
import spectral as sp

from src.group_evaluator import GroupEvaluator


# ============================================================================
# Paths
# ============================================================================

scratch = (
    r"C:\Users\Hyperspectral\Documents\GitHub"
    r"\python-tetracorder\scratch"
)

cube_hdr = scratch + r"\cuprite.95.cal.rtgc.v.hdr"
cube_file = scratch + r"\cuprite.95.cal.rtgc.v"

av95_hdr = scratch + r"\s07_AV95_envi.hdr"
av95_file = scratch + r"\s07_AV95_envi.sli"


# ============================================================================
# 1. Load Cuprite95 cube
# ============================================================================

data = sp.io.envi.open(
    cube_hdr,
    cube_file,
)

raw_cube = np.array(
    data.load(),
    dtype=np.float32,
)

metadata = sp.io.envi.read_envi_header(
    cube_hdr,
)

test_wvls = np.array(
    [float(x) for x in metadata["wavelength"]],
    dtype=float,
)

# Tetracorder Cuprite scale:
# DN 20000 = reflectance 1.0
test_cube = raw_cube * 0.00005

# Deleted-data sentinel, if present
test_cube[raw_cube == -32767] = np.nan


# ============================================================================
# 2. Load exact AVIRIS95.1 wavelength/FWHM records
# ============================================================================

av95 = sp.io.envi.open(
    av95_hdr,
    av95_file,
)

# SpectralLibrary.spectra shape:
# (number_of_spectra, number_of_bands)
av95_spectra = np.asarray(
    av95.spectra,
    dtype=float,
)

print("AV95 library shape:", av95_spectra.shape)
print("AV95 first records:", av95.names[:5])

# Record 0:
# Wavelengths in microns 224ch AVIRIS95.1
av95_wvls = av95_spectra[0]

# Record 1:
# Resolution in microns 224ch AVIRIS95.1
test_fwhm = av95_spectra[1]


# ============================================================================
# 3. Tetracorder AVIRIS95 deleted-channel mask
#
# Original 1-based channels:
#
# 1 2 13 31t33 81t83 95t97 108t111
# 153t167 172t175 224
# ============================================================================

valid_bands = np.ones(
    224,
    dtype=bool,
)

deleted_ranges = [
    (1, 2),
    (13, 13),
    (31, 33),
    (81, 83),
    (95, 97),
    (108, 111),
    (153, 167),
    (172, 175),
    (224, 224),
]

for start, stop in deleted_ranges:
    valid_bands[start - 1:stop] = False


# ============================================================================
# 4. Sanity checks
# ============================================================================

print("\nCube shape:", test_cube.shape)

print(
    "Cube reflectance range:",
    np.nanmin(test_cube),
    np.nanmean(test_cube),
    np.nanmax(test_cube),
)

print(
    "Cube wavelength range:",
    test_wvls.min(),
    test_wvls.max(),
)

print(
    "AV95 wavelength range:",
    av95_wvls.min(),
    av95_wvls.max(),
)

print(
    "FWHM range:",
    test_fwhm.min(),
    test_fwhm.max(),
)

print(
    "Bands:",
    len(test_wvls),
)

print(
    "Valid bands:",
    valid_bands.sum(),
)

print(
    "Deleted bands:",
    (~valid_bands).sum(),
)

delta = test_wvls - av95_wvls

print(
    "Max wavelength difference:",
    np.max(np.abs(delta)),
)

print(
    "Max VALID wavelength difference:",
    np.max(
        np.abs(
            delta[valid_bands]
        )
    ),
)

print(
    "Full wavelength order increasing:",
    np.all(
        np.diff(test_wvls) > 0
    ),
)

print(
    "Valid wavelength order increasing:",
    np.all(
        np.diff(
            test_wvls[valid_bands]
        ) > 0
    ),
)

assert test_cube.shape[-1] == 224
assert test_wvls.shape == (224,)
assert av95_wvls.shape == (224,)
assert test_fwhm.shape == (224,)
assert valid_bands.shape == (224,)


# ============================================================================
# 5. Small crop test
# ============================================================================

test_crop = test_cube[
    400:450,
    250:300,
    :
]

tetra = GroupEvaluator(
    test_cube,
    test_wvls,
    test_fwhm,
    target_valid_bands=valid_bands,
)

assert np.array_equal(
    tetra.evaluator.valid_bands,
    valid_bands,
)

print(
    "\nGroupEvaluator constructed successfully."
)


# ============================================================================
# 6. Inspect group winners
# ============================================================================

print("\nGroups:")

for gid, group in tetra.group_winners.items():

    if group is None:
        print(
            "group",
            gid,
            "-> None",
        )
        continue

    detected = (
        group["winner"] != None
    )

    print(
        "group",
        gid,
        "detections:",
        np.sum(detected),
        "max fit:",
        np.max(group["fit"]),
        "max depth:",
        np.max(group["depth"]),
        "max fd:",
        np.max(group["fit_depth"]),
    )


# ============================================================================
# 7. Inspect case winners
# ============================================================================

print("\nCases:")

for case_num, case in tetra.case_winners.items():

    if case is None:
        print(
            "case",
            case_num,
            "-> None",
        )
        continue

    detected = (
        case["winner"] != None
    )

    print(
        "case",
        case_num,
        "detections:",
        np.sum(detected),
        "max fit:",
        np.max(case["fit"]),
        "max depth:",
        np.max(case["depth"]),
        "max fd:",
        np.max(case["fit_depth"]),
    )


# ============================================================================
# 8. Explicit NVRES smoke test
# ============================================================================

evaluator = tetra.evaluator

for mid in (
    "red.edge.shift.1",
    "red.edge.shift.2",
):

    result = evaluator.evaluate(
        mid,
        test_crop,
    )

    print(
        "\n",
        mid,
    )

    if result is None:
        print("None")
        continue

    print(
        "max fit:",
        np.max(result["fit"]),
    )

    print(
        "max depth:",
        np.max(result["depth"]),
    )

    print(
        "max fd:",
        np.max(result["fit_depth"]),
    )

    print(
        "detections:",
        np.sum(
            result["fit"] > 0
        ),
    )
#%%
for gid, group in tetra.group_winners.items():

    if group is None:
        continue

    winners = group["winner"]

    names, counts = np.unique(
        winners[winners != None],
        return_counts=True,
    )

    if len(names) == 0:
        continue

    order = np.argsort(counts)[::-1]

    print(f"\nGROUP {gid}")

    for i in order[:10]:
        print(
            names[i],
            counts[i],
        )
#%%
import matplotlib.pyplot as plt
#tetra = GroupEvaluator(test_cube, test_wvls, test_fwhm)
group_winners = tetra.group_winners
for gid, group in group_winners.items():
    if group is not None:
        if np.any(group["fit"] > 0):
            plt.figure()
            plt.subplot(131)
            plt.title("fit")
            plt.imshow(group["fit"])
            
            plt.subplot(132)
            plt.title("depth")
            plt.imshow(group["depth"])
            
            plt.subplot(133)
            plt.title("fit_depth")
            plt.imshow(group["fit_depth"])
            plt.suptitle(f"Tetracorder group winners for: group {gid}")
            plt.savefig(f"C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/scratch/Tetracorder-group-winners-for-group-{gid}.png")
            plt.close()
    print(group.keys()) if group is not None else print(None)