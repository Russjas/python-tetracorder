
from pathlib import Path
import sys
import re
import gzip
import time
import struct
import importlib

import numpy as np

REPO = Path(r"C:\Users\Hyperspectral\Documents\GitHub\python-tetracorder") # editable
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
    
from tetracorderp.config import RRATIO_REFERENCES, NOGROUP0                      
from tetracorderp.material_evaluation import MaterialEvaluator          
from tetracorderp.group_evaluator import GroupEvaluator                 
#%% =============== directory and variable paths ==============================
# json rules parsed from native cmd files
RULES = REPO / "resources" / "tetracorder_rules_as_dict.json"
# Tetracorder preset mode
MODE = "default"
# python side reference database
REFERENCE_DB = REPO / "scratch" / "New_Convolver_work" / "tetracorder_rules_references_fwhm-out3.db"
# Root of the wsl holding the container holding the native Tetracorder
WSL_ROOT = Path(r"\\wsl.localhost\Ubuntu")
# setup and run directory of the baseline native run
NATIVE_RUN = WSL_ROOT / "home" / "hyperspectral" / "tetracorder-data" / "cuprite95" / "testrun1"
# Not all files used in the AVIRIS data runs on Tetracorder native are in the tetracorder-lite clone
# Some need to be accessed from the original repo instance
TET_ROOT = Path(r"C:\Users\Hyperspectral\Documents\GitHub\spectroscopy-tetracorder-extracted\spectroscopy-tetracorder-main")
# The SPECpr library the native run was configured with (restart r1-av95a).
# Supplies the sensor wavelength grid and bandpass widths for every test, and
# in the algorithmic test the convolved reference spectra themselves.
S06_AV95_TET =  TET_ROOT / "sl1" / "usgs" / "library06.conv" / "s06av95a"
# The convolved research library, for rules whose reference is [sprlb06].
R06_AV95_TET = TET_ROOT / "sl1" / "usgs" / "rlib06" / "r06av95a"
# Record 6 of s06av95a holds the AVIRIS95 wavelengths, per the native restart file.
WAVELENGTH_RECORD = 6
FWHM_RECORD = 12

# Image lines to evaluate. slice(None) is the whole scene; a narrow slice is
# for development only - the comparison reads native output at full size and
# slices it to match, so anything else gives a partial answer.
LINES = slice(None)
# This testrun calls the write_like_tetracorder() method, to produce
# native-formatted outputs. This directory specifies where they are written
PYTHON_OUTDIR = REPO / "scratch" / "test_runs" / "native_convolved"
#%% =============== helper functions to derive details from the native run ====
def wsl_path(path):
    """Map a Linux path recorded by Tetracorder native run onto the Windows WSL share."""
    path = str(path).strip()

    if path.startswith("/"):
        return WSL_ROOT / path.lstrip("/")

    return (NATIVE_RUN / path).resolve()

    
def read_envi_header(hdr_path):
    """Read the ENVI fields required to reconstruct the native cube."""

    text = Path(hdr_path).read_text(errors="replace")

    def get(key, default=None):
        match = re.search(
            rf"^\s*{key}\s*=\s*(.+)$",
            text,
            re.MULTILINE | re.IGNORECASE,
        )
        return match.group(1).strip() if match else default

    return {
        "samples": int(get("samples")),
        "lines": int(get("lines")),
        "bands": int(get("bands")),
        "offset": int(get("header offset", 0)),
        "data_type": int(get("data type")),
        "interleave": get("interleave").lower(),
        "byte_order": int(get("byte order", 0)),
    }


def read_native_cube(cube_path, dn_offset, deleted_dn, dn_scale):
    """
    Reproduce the native Tetracorder DN conversion:

        reflectance = float(DN + offset) * scale

    Deleted DN values become NaN.

    Returned shape is always:
        (lines, samples, bands)
    """

    hdr = read_envi_header(str(cube_path) + ".hdr")

    dtype = np.dtype({
        1: "u1",
        2: "i2",
        3: "i4",
        4: "f4",
        12: "u2",
    }[hdr["data_type"]])

    dtype = dtype.newbyteorder(
        ">" if hdr["byte_order"] == 1 else "<"
    )

    shape = {
        "bil": (hdr["lines"], hdr["bands"], hdr["samples"],),
        "bip": (hdr["lines"], hdr["samples"], hdr["bands"],),
        "bsq": (hdr["bands"], hdr["lines"], hdr["samples"],),}[hdr["interleave"]]

    raw = np.memmap(
        cube_path,
        dtype=dtype,
        mode="r",
        offset=hdr["offset"],
        shape=shape,
    )

    # Convert whatever the source interleave is to:
    # (lines, samples, bands)
    to_bip = {"bil": (0, 2, 1),
              "bip": (0, 1, 2),
              "bsq": (1, 2, 0),
              }[hdr["interleave"]]

    dn = np.transpose(raw, to_bip)

    # Do this explicitly in float32 to reproduce the native preparation
    cube = (dn.astype(np.float32) + np.float32(dn_offset)) * np.float32(dn_scale)

    cube[dn == deleted_dn] = np.nan

    return cube

def get_native_valid_bands(nchans):
    """
    Parse the native Tetracorder [DELETPTS] definition from cmds.start*
    and return a boolean mask:

        True  = usable band
        False = deleted band
    """

    deletpts = None

    for start_file in sorted(NATIVE_RUN.glob("cmds.start*")):

        with open(start_file, "r", encoding="utf-8", errors="replace",) as f:

            for line in f:

                if not line.startswith("==[DELETPTS]"):
                    continue

                value = (line[len("==[DELETPTS]"):].split("\\#")[0].strip())

                # Ignore unexpanded template entries
                if "DDDD" not in value:
                    deletpts = value
                    break

        if deletpts is not None:
            break

    if deletpts is None:
        raise ValueError("No usable ==[DELETPTS] definition found in native cmds.start files")

    deleted = np.zeros(nchans, dtype=bool,)

    for token in deletpts.split():

        if token.lower() == "c":
            break

        if "t" in token:

            start, stop = token.split("t")

            # Native channel numbers are 1-based and inclusive
            deleted[int(start) - 1: int(stop)] = True

        else:

            deleted[int(token) - 1] = True

    return ~deleted


def get_native_force_disabled():
    """
    Parse explicitly force-disabled groups and cases from:

        DISABLE/force-disable.txt
    """

    path = NATIVE_RUN / "DISABLE" / "force-disable.txt"

    if not path.exists():
        raise FileNotFoundError(
            f"Native force-disable file not found: {path}"
        )

    disabled_groups = set()
    disabled_cases = set()

    with open(path,"r", encoding="utf-8", errors="replace",) as f:
        for line in f:
            match = re.match(
                r"\s*DISABLE\s+(grp|cse)\s+(\d+)",
                line,
                re.IGNORECASE)

            if not match:
                continue

            number = int(match.group(2))

            if match.group(1).lower() == "grp":
                disabled_groups.add(number)

            else:
                disabled_cases.add(number)

    return disabled_groups, disabled_cases


def get_native_history_parameters():
    """
    Parse cube preparation parameters and scene physical conditions
    from the native Tetracorder history file.
    """

    with open(NATIVE_RUN / "history", "r", encoding="utf-8", errors="replace",) as f:
        history_lines = f.read().splitlines()

    result = {
        "cube_path": None,
        "dn_offset": None,
        "deleted_dn": None,
        "dn_scale": None,
        "temperature": None,
        "pressure": None,
    }

    for i, line in enumerate(history_lines):

        text = line.strip()

        # Cube path and DN conversion parameters
        if text.startswith("cube:"):

            result["cube_path"] = wsl_path(text[len("cube:"):].split()[0])

            tokens = history_lines[i + 1].split()

            result["dn_offset"] = int(tokens[0])
            result["deleted_dn"] = int(tokens[1])
            result["dn_scale"] = float(tokens[2])

            continue

        # Scene temperature
        match = re.match(
            r"temperature:?\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([CKck])\b",
            text,
        )

        if match:

            offset = (273.15 if match.group(3).lower() == "c" else 0.0)

            result["temperature"] = (
                float(match.group(1)) + offset,
                float(match.group(2)) + offset,
            )

            continue

        # Scene pressure
        match = re.match(
            r"pressure:?\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(\S+)",
            text,
        )

        if match:

            if not match.group(3).lower().startswith("bar"):
                raise ValueError(
                    f"Unexpected pressure unit: {match.group(3)}"
                )

            result["pressure"] = (
                float(match.group(1)),
                float(match.group(2)),
            )

    missing = [key for key, value in result.items() if value is None]

    if missing:
        raise ValueError("Missing native history values: " + ", ".join(missing))

    return result
# =============== SPECpr reading =============================================

SPEC_PR_RECORD_BYTES = 1536

DATA_HEADER_FORMAT = ">i40s8s16i60s74s74s74s74s6i260f"
DATA_CONTINUATION_FORMAT = ">i383f"


def _read_physical_record(path, record_number):
    with open(path, "rb") as f:
        f.seek(record_number * SPEC_PR_RECORD_BYTES)
        record = f.read(SPEC_PR_RECORD_BYTES)

    if len(record) != SPEC_PR_RECORD_BYTES:
        raise EOFError(
            f"Could not read SPECpr record {record_number} from {path}"
        )

    return record


def read_specpr(path, record_number):
    """Read one SPECpr data record, including continuation records."""

    first = struct.unpack(DATA_HEADER_FORMAT, _read_physical_record(path, record_number),)

    if first[0] % 4 != 0:
        raise ValueError(f"{Path(path).name} record {record_number} is not a data record")

    title = first[1].decode("latin-1").rstrip("\x00 ")

    nchans = first[10]

    data = np.empty(nchans, dtype=np.float32,)

    written = min(nchans, 256,)

    data[:written] = first[34:34 + written]

    record_number += 1

    while written < nchans:

        continuation = struct.unpack(
            DATA_CONTINUATION_FORMAT,
            _read_physical_record(path, record_number),
        )

        if continuation[0] % 4 != 1:
            raise ValueError(f"{Path(path).name}: expected continuation record at {record_number}")

        count = min(383, nchans - written,)

        data[written:written + count] = continuation[1:1 + count]

        written += count
        record_number += 1

    # SPECpr deleted-point sentinel
    data[data < -1.0e34] = np.nan

    return title, data

#%% =============== Find the cube and parameters used in the native run =======
# parse the history of the native run to reuse the cube and parameters
native = get_native_history_parameters()
cube_path = native["cube_path"]
dn_offset = native["dn_offset"]
deleted_dn = native["deleted_dn"]
dn_scale = native["dn_scale"]
# produce an ndarray prepared in the same way as native Tetracorder preprocesses
prepared_test_cube = read_native_cube(cube_path, dn_offset, deleted_dn, dn_scale)

# read the sensor specific wavelengths and fwhm, stored in libraries distributed 
# with spectroscopy-tetracorder
wavelength_title, wavelengths = read_specpr(S06_AV95_TET, WAVELENGTH_RECORD,)
fwhm_title, fwhm = read_specpr(S06_AV95_TET, FWHM_RECORD,)
# construct the bad band list/valid bands used for the cube in the native run
valid_bands = get_native_valid_bands(wavelengths.size)

# Native Tetracorder has scence physical constraint parameter we need to get
# they can control material resolution behaviour
scene_temperature = native["temperature"]
scene_pressure = native["pressure"]

# The orinal Cuprite95 build benchmark on spectroscopy-tetracorder had 
# force disabled some groups and cases
disabled_groups, disabled_cases = get_native_force_disabled()

print("wavelength record title:", wavelength_title)
print("wavelength record shape:", wavelengths.shape)
print("fwhm record title:", fwhm_title)
print("fwhm record shape:", fwhm.shape)
print("cube:", prepared_test_cube.shape)
print("dtype:", prepared_test_cube.dtype)
print("cube shape:", prepared_test_cube.shape)
print("wavelengths/fwhm/cube bands match:", 
      wavelengths.shape == fwhm.shape == (prepared_test_cube.shape[-1],))


#%% ================ perform the python run ===================================

t0 = time.perf_counter()
run = GroupEvaluator(
    prepared_test_cube, wavelengths, fwhm, mode=MODE,
    target_valid_bands=valid_bands,
    reference_file=str(REFERENCE_DB), rules_file=str(RULES),
    disabled_groups=disabled_groups, disabled_cases=disabled_cases,
    temperature=scene_temperature, pressure=scene_pressure,
)
print(f"\npackage run: {time.perf_counter() - t0:.0f} s")
print(f"materials disabled: {len(run.disabled_materials)}")
for mid in sorted(run.disabled_materials):
    rule = run.evaluator.rules[mid]
    print(f"    {mid}  ({rule['kind']} {rule['number']})")
rules = run.evaluator.rules


#%% =============helper functions for native output comparison ================

def get_output_dirs(root):

    return {
        path.name: path
        for pattern in ("group.*", "case.*")
        for path in root.glob(pattern)
        if path.is_dir()
    }

def read_vicar_image(path):

    with gzip.open(path, "rb") as f:
        raw = f.read()

    match = re.match(rb"LBLSIZE=(\d+)", raw)

    if match is None:
        raise ValueError(f"No VICAR LBLSIZE in {path}")

    lblsiz = int(match.group(1))

    return np.frombuffer(
        raw[lblsiz:],
        dtype=np.uint8,
    )

def load_native_style_output(directory):

    directory = Path(directory)

    outputs = {}

    for fit_file in directory.glob("*.fit.gz"):

        basename = fit_file.name.removesuffix(".fit.gz")

        depth_file = directory / f"{basename}.depth.gz"
        fd_file = directory / f"{basename}.fd.gz"

        outputs[basename] = {
            "fit": read_vicar_image(fit_file),
            "depth": read_vicar_image(depth_file),
            "fd": read_vicar_image(fd_file),
        }

    return outputs
#========== helpers for statistic calculation =================================

def stats(native, python, mask):
    """
    Compare two native-style uint8 outputs over a supplied boolean mask.

    The mask defines the pixels considered valid for comparison in both
    outputs, e.g. their common classified pixels, optionally further
    restricted to nonzero depth or fit-depth values.
    """

    delta = (python[mask].astype(np.int16) - native[mask].astype(np.int16))

    if delta.size == 0:
        return {
            "n": 0,
            "exact": np.nan,
            "within1": np.nan,
            "median": np.nan,
            "max": 0,
        }

    return {"n": delta.size,"exact": np.mean(delta == 0),
            "within1": np.mean(np.abs(delta) <= 1),
            "median": float(np.median(delta)),
            "max": int(np.abs(delta).max()),
            }

def compare_material(material, py, nat):

    nat_present = nat["fit"] > 0
    py_present = py["fit"] > 0

    both = nat_present & py_present
    union = nat_present | py_present

    return {
        "material": material,

        "native": int(nat_present.sum()),
        "python": int(py_present.sum()),
        "both": int(both.sum()),
        "union": int(union.sum()),

        "python_only": int((py_present & ~nat_present).sum()),

        "native_only": int((nat_present & ~py_present).sum()),

        "jaccard": (both.sum() / union.sum() if union.any() else np.nan),

        "fit": stats(nat["fit"], py["fit"], both,),

        "depth": stats(nat["depth"], py["depth"], 
                       both & (nat["depth"] > 0) & (py["depth"] > 0)),

        "fd": stats(nat["fd"], py["fd"],
            both & (nat["fd"] > 0) & (py["fd"] > 0)),
    }

def compare_group(py, nat):

    rows = []

    common_materials = sorted(set(py) & set(nat))

    for material in common_materials:
        rows.append(compare_material(
                material,
                py[material],
                nat[material]))

    return rows

def summarize_group(rows):

    both = sum(row["both"] for row in rows)

    union = sum(row["union"] for row in rows)

    fit_n = sum(row["fit"]["n"] for row in rows)

    fit_exact = sum(row["fit"]["exact"] * row["fit"]["n"]
        for row in rows if row["fit"]["n"])

    native = sum(row["native"] for row in rows)

    python = sum(row["python"] for row in rows)

    return (
        both / union if union else np.nan,
        fit_exact / fit_n if fit_n else np.nan,
        native,
        python,
    )

def format_report(group_results, summary):
    """
    Format per-group comparison results and the final summary as Markdown.

    Returns the complete report as a single string suitable for both
    console output and writing to a .md file.
    """

    lines = []

    for dirname, rows in group_results.items():

        lines.append(f"## {dirname}")
        lines.append("")

        lines.append(
            "| Material | Native | Python | Jaccard | "
            "Fit exact | Depth exact | Fit max |"
        )
        lines.append(
            "|---|---:|---:|---:|---:|---:|---:|"
        )

        for row in sorted(
            rows,
            key=lambda row: -(row["native"] + row["python"]),
        ):

            if row["native"] == 0 and row["python"] == 0:
                continue

            lines.append(
                f"| `{row['material']}` "
                f"| {row['native']:,} "
                f"| {row['python']:,} "
                f"| {row['jaccard']:.3f} "
                f"| {row['fit']['exact']:.3f} "
                f"| {row['depth']['exact']:.3f} "
                f"| {row['fit']['max']} |"
            )

        lines.append("")

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------

    lines.append("## Summary")
    lines.append("")

    lines.append(
        "| Group | Materials | Native only | Python only | "
        "Native pixels | Python pixels | Jaccard | Fit exact |"
    )
    lines.append(
        "|---|---:|---:|---:|---:|---:|---:|---:|"
    )

    for row in summary:

        lines.append(
            f"| {row['group']} "
            f"| {row['materials']} "
            f"| {row['native_only']} "
            f"| {row['python_only']} "
            f"| {row['native_pixels']:,} "
            f"| {row['python_pixels']:,} "
            f"| {row['jaccard']:.4f} "
            f"| {row['fit_exact']:.4f} |"
        )

    return "\n".join(lines)

#%% =====--- compare python and native outputs ================================
# write the tetracorder style outputs from the completed python run
run.write_like_tetracorder(PYTHON_OUTDIR)

native_output_dirs = get_output_dirs(NATIVE_RUN)
python_output_dirs = get_output_dirs(PYTHON_OUTDIR)

print("native output dirs:", len(native_output_dirs))
print("python output dirs:", len(python_output_dirs))

native_only_dirs = (set(native_output_dirs) - set(python_output_dirs))

python_only_dirs = (set(python_output_dirs) - set(native_output_dirs))

common_dirs = (set(native_output_dirs) & set(python_output_dirs))

print("common dirs:", len(common_dirs))
print("native only:", sorted(native_only_dirs))
print("python only:", sorted(python_only_dirs))

group_results = {}
summary = []

for dirname in sorted(common_dirs):

    nat = load_native_style_output(native_output_dirs[dirname])

    py = load_native_style_output(python_output_dirs[dirname])

    rows = compare_group(py, nat)

    group_results[dirname] = rows

    jaccard, fit_exact, native_pixels, python_pixels = summarize_group(rows)

    native_only = set(nat) - set(py)
    python_only = set(py) - set(nat)

    summary.append({
        "group": dirname,
        "materials": len(rows),
        "native_only": len(native_only),
        "python_only": len(python_only),
        "native_pixels": native_pixels,
        "python_pixels": python_pixels,
        "jaccard": jaccard,
        "fit_exact": fit_exact,
    })
report = format_report(group_results, summary)

# Console
print(report)

# Markdown
REPORT_FILE = PYTHON_OUTDIR / "fidelity_test3.md"

with open(REPORT_FILE, "w", encoding="utf-8",) as f:
    f.write(report)

print(f"\nReport written to: {REPORT_FILE}")

#%% =========== End of test, anything below here is messing =====================
