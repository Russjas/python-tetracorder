"""
Run the python-tetracorder package end to end on the native cuprite95 test
inputs and compare with the native Tetracorder 6.00 run.

This script only PREPARES inputs:
    cube          native int16 cube, (DN + offset) * scale, deleted DN -> NaN
    wavelengths   native wavelength record (s06av95a record 6)
    valid bands   native [DELETPTS]


and hands them to the package. Everything else is the package itself:
GroupEvaluator.evaluate() (groups, group 0, NOGROUP0, cases) and
MaterialEvaluator (features, tests, weights, constraints, NOTs).
The references are from the SQLite db shipped with the package, extracted from the 
preconvolution specpr Tetracorder libraries, and then convolved using the vendored Tetrapy 
convolver

"""
from pathlib import Path
import sys
import re
import gzip
import time
import struct
import importlib

import numpy as np


# =============================================================================
# PATHS  (edit these)
# =============================================================================

REPO = Path(r"C:\Users\Hyperspectral\Documents\GitHub\python-tetracorder")
RULES = REPO / "resources" / "tetracorder_rules_as_dict.json"
MODE = "default"
REFERENCE_DB = REPO / "scratch" / "New_Convolver_work" / "tetracorder_rules_references_fwhm-out3.db"

WSL_ROOT = Path(r"\\wsl.localhost\Ubuntu")
NATIVE_RUN = WSL_ROOT / "home" / "hyperspectral" / "tetracorder-data" / "cuprite95" / "testrun1"

TET_ROOT = Path(
    r"C:\Users\Hyperspectral\Documents\GitHub"
    r"\spectroscopy-tetracorder-extracted"
    r"\spectroscopy-tetracorder-main"
)

# Groups / cases to compare against native output ("all" or a tuple of numbers).
# All groups and cases are always evaluated, as GroupEvaluator does.
# Group 0 materials are compared inside every group that includes group 0,
# as native writes a group-0 output file into each such group directory.
COMPARE_GROUPS = "all"
COMPARE_CASES = "all"
DETAIL = ("calcite.ws272.g2", "aragonite")
SHOW_MATERIALS = True        # per-material table under each group / case

# Evaluate a block of image lines, e.g. slice(0, 100); slice(None) = all.
LINES = slice(None)

CUBE_CACHE_DIR = REPO / "scratch" / "native_cube_cache"
WAVELENGTH_RECORD = 6                      # restart r1-av95a: itrol = Y 6

# Fallbacks if the native run files cannot be read
DEFAULT_DELETPTS = "1 2 13 31t33 81t83 95t97 108t111 153t167 172t175 224 c"
DEFAULT_OFFSET, DEFAULT_DELETED_DN, DEFAULT_SCALE = 0, -32767, 0.00005
DEFAULT_DEPTH_FULL_SCALE = 0.5


def first_existing(*candidates):
    for candidate in candidates:
        if Path(candidate).is_file():
            return Path(candidate)
    raise FileNotFoundError("none of these exist:\n  " +
                            "\n  ".join(str(c) for c in candidates))


# restart r1-av95a:  y = /sl1/usgs/library06.conv/s06av95a
#                    w = /sl1/usgs/rlib06/r06av95a
S06_AV95 = first_existing(
    WSL_ROOT / "sl1" / "usgs" / "library06.conv" / "s06av95a",
    TET_ROOT / "sl1" / "usgs" / "library06.conv" / "s06av95a",
)

print("s06av95a:", S06_AV95)



# =============================================================================
# PACKAGE  (reloaded so edits to tetracorderp are picked up in the same kernel)
# =============================================================================

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import tetracorderp.config                                              # noqa: E402
import tetracorderp.tetracorder_ops                                     # noqa: E402
import tetracorderp.material_evaluation                                 # noqa: E402
import tetracorderp.group_evaluator                                     # noqa: E402

for module in (tetracorderp.config, tetracorderp.tetracorder_ops,
               tetracorderp.material_evaluation, tetracorderp.group_evaluator):
    importlib.reload(module)

from tetracorderp.config import RRATIO_REFERENCES                       # noqa: E402
from tetracorderp.material_evaluation import MaterialEvaluator          # noqa: E402
from tetracorderp.group_evaluator import GroupEvaluator                 # noqa: E402


# =============================================================================
# SPECpr reader
# =============================================================================

RECORD_BYTES = 1536
HEADER_FORMAT = ">i40s8s16i60s74s74s74s74s6i260f"
CONTINUATION_FORMAT = ">i383f"
_LIBRARY_BYTES = {}


def _physical_record(path, record_number):
    path = Path(path)
    if path not in _LIBRARY_BYTES:
        _LIBRARY_BYTES[path] = path.read_bytes()
    start = record_number * RECORD_BYTES
    record = _LIBRARY_BYTES[path][start:start + RECORD_BYTES]
    if len(record) != RECORD_BYTES:
        raise EOFError(f"Could not read SPECpr record {record_number} from {path}")
    return record


def read_specpr(path, record_number):
    """One SPECpr data record incl. continuations; deleted points -> NaN."""
    record_number = int(record_number)
    first = struct.unpack(HEADER_FORMAT, _physical_record(path, record_number))
    if first[0] % 4 != 0:
        raise ValueError(f"{Path(path).name} rec {record_number}: not a data record")

    title = first[1].decode("latin-1").rstrip("\x00 ")
    nch = first[10]
    data = np.empty(nch, dtype=np.float32)
    written = min(nch, 256)
    data[:written] = first[34:34 + written]
    rec = record_number + 1
    while written < nch:
        cont = struct.unpack(CONTINUATION_FORMAT, _physical_record(path, rec))
        if cont[0] % 4 != 1:
            raise ValueError(f"{Path(path).name}: expected continuation at {rec}")
        count = min(383, nch - written)
        data[written:written + count] = cont[1:1 + count]
        written += count
        rec += 1

    data[data < -1.0e34] = np.nan                  # SPECpr deleted point
    return title, data


# =============================================================================
# NATIVE RUN SETTINGS
# =============================================================================

def native_start_lines():
    for start in sorted(NATIVE_RUN.glob("cmds.start*")):
        yield from start.read_text(errors="replace").splitlines()


def native_deletpts():
    for line in native_start_lines():
        if line.startswith("==[DELETPTS]"):
            value = line[len("==[DELETPTS]"):].split("\\#")[0].strip()
            if "DDDD" not in value:
                return value
    return DEFAULT_DELETPTS


def native_output_dirs():
    """==[DIRg2]group.2um/ -> {("group", 2): ...}; ==[DIRc6]case.x/ -> {("case", 6): ...}"""
    dirs = {}
    for line in native_start_lines():
        m = re.match(r"==\[DIR([gc])(\d+)\](\S+)", line)
        if m:
            kind = "group" if m.group(1) == "g" else "case"
            dirs[(kind, int(m.group(2)))] = NATIVE_RUN / m.group(3).strip("/")
    return dirs


def parse_deletpts(text, nchans):
    deleted = np.zeros(nchans, dtype=bool)
    for token in text.split():
        if token.lower() == "c":
            break
        if "t" in token:
            a, b = token.split("t")
            deleted[int(a) - 1:int(b)] = True
        else:
            deleted[int(token) - 1] = True
    return ~deleted


def native_cube_parameters():
    history = NATIVE_RUN / "history"
    if history.exists():
        lines = history.read_text(errors="replace").splitlines()
        for i, line in enumerate(lines):
            if line.startswith("cube:"):
                path = line[len("cube:"):].split()[0]
                path = (WSL_ROOT / path.lstrip("/")) if path.startswith("/") \
                    else (NATIVE_RUN / path).resolve()
                t = lines[i + 1].split()
                return path, int(t[0]), int(t[1]), float(t[2])
    return (NATIVE_RUN.parent / "cube" / "cuprite.95.cal.rtgc.v",
            DEFAULT_OFFSET, DEFAULT_DELETED_DN, DEFAULT_SCALE)


def native_depth_scales():
    scales = {}
    path = NATIVE_RUN / "AAA.info" / "material-DN-scalling.txt"
    if path.exists():
        for line in path.read_text(errors="replace").splitlines():
            m = re.search(r"(group|case)\s+(-?\d+)\s+(\S+)\s+DN=\s*(\d+)\s*=\s*([0-9.]+)", line)
            if m:
                scales[(m.group(1), int(m.group(2)), m.group(3))] = (
                    int(m.group(4)) / float(m.group(5)))
    return scales


def native_conditions():
    """(temperature K, pressure bar) of the data set, from cmds.start."""
    temperature = pressure = None
    for line in native_start_lines():
        text = line.split("\\#")[0].strip()
        m = re.match(r"temperature:?\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([CKck])\b", text)
        if m:
            offset = 273.15 if m.group(3) in "Cc" else 0.0
            temperature = (float(m.group(1)) + offset, float(m.group(2)) + offset)
        m = re.match(r"pressure:?\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(\S+)", text)
        if m:
            if not m.group(3).lower().startswith("bar"):
                raise ValueError(f"unexpected pressure unit: {m.group(3)}")
            pressure = (float(m.group(1)), float(m.group(2)))
    return temperature, pressure


def native_force_disabled():
    """(groups, cases) from DISABLE/force-disable.txt."""
    groups, cases = set(), set()
    path = NATIVE_RUN / "DISABLE" / "force-disable.txt"
    if path.exists():
        for line in path.read_text(errors="replace").splitlines():
            m = re.match(r"\s*DISABLE\s+(grp|cse)\s+(\d+)", line.strip(), re.IGNORECASE)
            if m:
                (groups if m.group(1).lower() == "grp" else cases).add(int(m.group(2)))
    return groups, cases

def read_envi_header(hdr_path):
    text = Path(hdr_path).read_text(errors="replace")

    def get(key, default=None):
        m = re.search(rf"^\s*{key}\s*=\s*(.+)$", text, re.MULTILINE | re.IGNORECASE)
        return m.group(1).strip() if m else default

    return dict(samples=int(get("samples")), lines=int(get("lines")),
                bands=int(get("bands")), offset=int(get("header offset", 0)),
                data_type=int(get("data type")), interleave=get("interleave").lower(),
                byte_order=int(get("byte order", 0)))


def read_native_cube(cube_path, offset, deleted_dn, scale):
    """read2sheet.r: float(DN + offset) * scale; DN == deleted -> deleted."""
    hdr = read_envi_header(str(cube_path) + ".hdr")
    dtype = np.dtype({1: "u1", 2: "i2", 3: "i4", 4: "f4", 12: "u2"}[hdr["data_type"]])
    dtype = dtype.newbyteorder(">" if hdr["byte_order"] == 1 else "<")
    shape = {"bil": (hdr["lines"], hdr["bands"], hdr["samples"]),
             "bip": (hdr["lines"], hdr["samples"], hdr["bands"]),
             "bsq": (hdr["bands"], hdr["lines"], hdr["samples"])}[hdr["interleave"]]
    raw = np.memmap(cube_path, dtype=dtype, mode="r", offset=hdr["offset"], shape=shape)
    dn = np.transpose(raw, {"bil": (0, 2, 1), "bip": (0, 1, 2),
                            "bsq": (1, 2, 0)}[hdr["interleave"]])
    cube = (dn.astype(np.float32) + np.float32(offset)) * np.float32(scale)
    cube[dn == deleted_dn] = np.nan
    return cube


# =============================================================================
# PREPARE INPUTS
# =============================================================================



wl_title, wavelengths = read_specpr(S06_AV95, WAVELENGTH_RECORD)
print("wavelengths:", wl_title)

valid_bands = parse_deletpts(native_deletpts(), wavelengths.size)
print("valid channels:", int(valid_bands.sum()))

fwhm = None
for rec in range(1, 60):
    try:
        title, data = read_specpr(S06_AV95, rec)
    except (ValueError, EOFError, struct.error):
        continue
    if (title.startswith(("Resolution", "Bandpass")) or "FWHM" in title) \
            and data.size == wavelengths.size:
        fwhm = data
        print(f"sensor FWHM: rec {rec} {title!r}  "
              f"{np.nanmin(fwhm):.4f}-{np.nanmax(fwhm):.4f}")
        break
if fwhm is None:
    raise RuntimeError("no resolution record found in s06av95a")

cube_path, dn_offset, deleted_dn, dn_scale = native_cube_parameters()
print(f"cube: {cube_path}  offset={dn_offset} deleted={deleted_dn} scale={dn_scale}")

CUBE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
cache_file = CUBE_CACHE_DIR / (
    f"{Path(cube_path).name}_off{dn_offset}_del{deleted_dn}_sc{dn_scale:g}.npy")
if cache_file.exists():
    full_cube = np.load(cache_file)
else:
    full_cube = read_native_cube(cube_path, dn_offset, deleted_dn, dn_scale)
    np.save(cache_file, full_cube)

FULL_SHAPE = full_cube.shape[:2]
cube = np.ascontiguousarray(full_cube[LINES])
del full_cube
print("cube:", cube.shape, cube.dtype)


def native_rule_spectra(rules):
    spectra, failures = {}, []
    for rule_id, rule in rules.items():
        small = rule["library_records"]["SMALL"]
        try:
            _, data = read_specpr(LIBRARY_FILES[small["library"]], small["record"])
        except Exception as exc:
            failures.append(f"{rule_id}: {small} {exc!r}")
            continue
        if data.shape != wavelengths.shape:
            failures.append(f"{rule_id}: {data.shape} channels")
            continue
        spectra[rule_id] = data
    if failures:
        raise RuntimeError("native references not loaded:\n  " + "\n  ".join(failures))
    return spectra


#%% =============================================================================
# RUN THE PACKAGE
# =============================================================================
TEMPERATURE, PRESSURE = native_conditions()
FORCED_GROUPS, FORCED_CASES = native_force_disabled()
print(f"data conditions: {TEMPERATURE} K, {PRESSURE} bar")
print(f"force-disabled: groups {sorted(FORCED_GROUPS)}, cases {sorted(FORCED_CASES)}")

t0 = time.perf_counter()
run = GroupEvaluator(
    cube, wavelengths, fwhm, mode=MODE,
    target_valid_bands=valid_bands,
    reference_file=str(REFERENCE_DB), rules_file=str(RULES),
    disabled_groups=FORCED_GROUPS, disabled_cases=FORCED_CASES,
    temperature=TEMPERATURE, pressure=PRESSURE,
)
print(f"\npackage run: {time.perf_counter() - t0:.0f} s")
print(f"materials disabled: {len(run.disabled_materials)}")
for mid in sorted(run.disabled_materials):
    rule = run.evaluator.rules[mid]
    print(f"    {mid}  ({rule['kind']} {rule['number']})")
rules = run.evaluator.rules


#% =============================================================================
# COMPARE WITH NATIVE OUTPUT
# =============================================================================
import importlib, tetracorderp.tetracorder_ops, tetracorderp.material_evaluation
importlib.reload(tetracorderp.tetracorder_ops)
importlib.reload(tetracorderp.material_evaluation)
from tetracorderp.config import NOGROUP0                                # noqa: E402


def nint_uint8(values, scale):
    """cublineout.r: nint(x * scale), <= 0 -> 0, > 255 -> 255."""
    x = np.asarray(np.ma.filled(np.ma.asarray(values, dtype=np.float64), 0.0)) * scale
    x = np.where(np.isfinite(x), x, 0.0)
    return np.clip(np.floor(x + 0.5), 0, 255).astype(np.uint8)


def read_native_byte_gz(path):
    n = int(np.prod(FULL_SHAPE))
    with gzip.open(path, "rb") as f:
        raw = f.read()
    return np.frombuffer(raw[len(raw) - n:], dtype=np.uint8).reshape(FULL_SHAPE)[LINES]


def output_info(rule):
    """(basename, bits) from the rule's output block."""
    lines = [l.split("\\#")[0].strip() for l in rule.get("output_raw", "").splitlines()]
    lines = [l for l in lines if l]
    if len(lines) < 3:
        return None, None
    return lines[1].split()[0], int(lines[2].split()[0])


def stats(native, python, mask):
    d = python[mask].astype(np.int16) - native[mask].astype(np.int16)
    if d.size == 0:
        return dict(n=0, exact=np.nan, within1=np.nan, median=np.nan, max=0)
    return dict(n=d.size, exact=np.mean(d == 0), within1=np.mean(np.abs(d) <= 1),
                median=float(np.median(d)), max=int(np.abs(d).max()))


output_dirs = native_output_dirs()
depth_scales = native_depth_scales()


def compare_set(kind, number, resolved, members):
    """Compare every member material of one group / case with native."""
    out_dir = output_dirs.get((kind, number))
    rows, missing = [], 0
    shape = cube.shape[:2]
    winner = resolved["winner"] if resolved is not None else np.full(shape, None, object)
    zeros = np.zeros(shape)

    for mid in members:
        rule = rules[mid]
        basename, bits = output_info(rule)
        if out_dir is None or basename is None or bits != 8:
            missing += 1
            continue
        files = [out_dir / f"{basename}.{k}.gz" for k in ("fit", "depth", "fd")]
        if not all(f.exists() for f in files):
            missing += 1
            continue

        n_fit, n_depth, n_fd = (read_native_byte_gz(f) for f in files)
        depth_scale = depth_scales.get(
            (rule["kind"], int(rule["number"]), basename),
            255 / DEFAULT_DEPTH_FULL_SCALE)

        mine = winner == mid
        get = (lambda k: resolved[k]) if resolved is not None else (lambda k: zeros)
        p_fit = nint_uint8(np.where(mine, get("fit"), 0.0), 255.0)
        p_depth = nint_uint8(np.where(mine, get("depth"), 0.0), depth_scale)
        p_fd = nint_uint8(np.where(mine, get("fit_depth"), 0.0), depth_scale)

        nat, pyt = n_fit > 0, p_fit > 0
        both = nat & pyt
        union = nat | pyt
        row = dict(
            material=mid, group0=(rule["kind"] == "group" and int(rule["number"]) == 0),
            native=int(nat.sum()), python=int(pyt.sum()), both=int(both.sum()),
            union=int(union.sum()),
            jaccard=(both.sum() / union.sum()) if union.any() else np.nan,
            fit=stats(n_fit, p_fit, both),
            depth=stats(n_depth, p_depth, both & (n_depth > 0) & (p_depth > 0)),
        )
        rows.append(row)

        if mid in DETAIL:
            print(f"\n  --- {mid} ({kind} {number}, {basename})")
            print(f"      native {row['native']}  python {row['python']}  both {row['both']}"
                  f"  jaccard {row['jaccard']:.4f}")
            no = nat & ~pyt
            if no.any():
                ids, counts = np.unique(winner[no].astype(str), return_counts=True)
                top = np.argsort(counts)[::-1][:5]
                print("      Python winner on native-only pixels:",
                      {str(ids[i]): int(counts[i]) for i in top})
            for k in ("fit", "depth"):
                st = row[k]
                print(f"      common {k:5s} n={st['n']} exact={st['exact']:.4f} "
                      f"max={st['max']}")
    return rows, missing


def print_table(rows):
    print(f"    {'material':42s} {'native':>7s} {'python':>7s} {'jacc':>6s} "
          f"{'fit=':>6s} {'dep=':>6s} {'fitmax':>6s}")
    for r in sorted(rows, key=lambda r: -(r["native"] + r["python"])):
        if r["native"] == 0 and r["python"] == 0:
            continue
        name = ("[g0] " if r["group0"] else "") + r["material"]
        print(f"    {name[:42]:42s} {r['native']:7d} {r['python']:7d} "
              f"{r['jaccard']:6.3f} {r['fit']['exact']:6.3f} "
              f"{r['depth']['exact']:6.3f} {r['fit']['max']:6d}")


def pooled(rows):
    both = sum(r["both"] for r in rows)
    union = sum(r["union"] for r in rows)
    exact_n = sum(r["fit"]["n"] for r in rows)
    exact = sum(r["fit"]["exact"] * r["fit"]["n"] for r in rows if r["fit"]["n"])
    return (both / union if union else np.nan,
            exact / exact_n if exact_n else np.nan,
            sum(r["native"] for r in rows), sum(r["python"] for r in rows))


group_numbers = sorted({int(r["number"]) for r in rules.values()
                        if r["kind"] == "group" and int(r["number"]) > 0})
case_numbers = sorted({int(r["number"]) for r in rules.values() if r["kind"] == "case"})

if COMPARE_GROUPS != "all":
    group_numbers = [g for g in group_numbers if g in COMPARE_GROUPS]
if COMPARE_CASES != "all":
    case_numbers = [c for c in case_numbers if c in COMPARE_CASES]
group_numbers = [g for g in group_numbers if g not in run.disabled_groups]
case_numbers = [c for c in case_numbers if c not in run.disabled_cases]

group0_members = [m for m, r in rules.items()
                  if r["kind"] == "group" and int(r["number"]) == 0]
summary = []

for g in group_numbers:
    members = [m for m, r in rules.items()
               if r["kind"] == "group" and int(r["number"]) == g]
    if g not in NOGROUP0:
        members += group0_members
    print("\n" + "=" * 78 + f"\nGROUP {g}  ({output_dirs.get(('group', g))})\n" + "=" * 78)
    rows, missing = compare_set("group", g, run.group_winners.get(g), members)
    jac, fexact, nn, pn = pooled(rows)
    summary.append((f"group {g}", len(rows), missing, nn, pn, jac, fexact))
    print(f"  compared {len(rows)} materials ({missing} without native output); "
          f"native px {nn}, python px {pn}, pooled jaccard {jac:.4f}, "
          f"common-pixel fit exact {fexact:.4f}")
    if SHOW_MATERIALS:
        print_table(rows)

for c in case_numbers:
    members = [m for m, r in rules.items()
               if r["kind"] == "case" and int(r["number"]) == c]
    print("\n" + "=" * 78 + f"\nCASE {c}  ({output_dirs.get(('case', c))})\n" + "=" * 78)
    rows, missing = compare_set("case", c, run.case_winners.get(c), members)
    jac, fexact, nn, pn = pooled(rows)
    summary.append((f"case {c}", len(rows), missing, nn, pn, jac, fexact))
    print(f"  compared {len(rows)} materials ({missing} without native output); "
          f"native px {nn}, python px {pn}, pooled jaccard {jac:.4f}, "
          f"common-pixel fit exact {fexact:.4f}")
    if SHOW_MATERIALS:
        print_table(rows)

print("\n" + "=" * 78 + "\nSUMMARY\n" + "=" * 78)
print(f"{'set':10s} {'mats':>5s} {'noout':>5s} {'native':>8s} {'python':>8s} "
      f"{'jacc':>7s} {'fit=':>7s}")
for label, n, miss, nn, pn, jac, fexact in summary:
    print(f"{label:10s} {n:5d} {miss:5d} {nn:8d} {pn:8d} {jac:7.4f} {fexact:7.4f}")
    
#%%
run.write_npz(r"C:\Users\Hyperspectral\Documents\GitHub\python-tetracorder\tests\Tetrapy_convolver_full_run_results\tet_conv_full.npz")
run.write_like_tetracorder(r"C:\Users\Hyperspectral\Documents\GitHub\python-tetracorder\tests\Tetrapy_convolver_full_run_results\tet_like_output")
run.write_jpgs(r"C:\Users\Hyperspectral\Documents\GitHub\python-tetracorder\tests\Tetrapy_convolver_full_run_results\jpgs_output")


#%%END TEST SCRIPT, MESSING BELOW
#%% profile a strip
import cProfile, pstats
from tetracorderp.group_evaluator import GroupEvaluator
LINES = slice(0, 100)
cube = np.ascontiguousarray(np.load(cache_file)[LINES])
cProfile.run("GroupEvaluator(cube, wavelengths, fwhm, mode=MODE, target_valid_bands=valid_bands, "
             "reference_file=str(REFERENCE_DB), rules_file=str(RULES), disabled_groups=FORCED_GROUPS, "
             "disabled_cases=FORCED_CASES, temperature=TEMPERATURE, pressure=PRESSURE)", "prof.out")
pstats.Stats("prof.out").sort_stats("cumulative").print_stats(25)
#%%
#%% lawn grass: old ascii copy vs splib06a copy, both against native 7644
#%% [RATIOGVEG1] lawn grass: splib06a copy vs old ASCII copy, both against native s06av95a 7644
import sqlite3
from tetracorderp.convolve import Convolver

_, native = read_specpr(S06_AV95, 7644)
conv = Convolver(wavelengths, fwhm)

def lawn_grass(db):
    con = sqlite3.connect(db)
    cols = [r[1] for r in con.execute("PRAGMA table_info(Spectra)")]
    sel = "XData, YData" + (", FData" if "FData" in cols else "")
    row = con.execute(f"SELECT {sel} FROM Samples JOIN Spectra USING (SampleID) "
                      "WHERE Library = 'splib06b' AND ConvolvedRecord = 7644").fetchone()
    con.close()
    arrays = [np.frombuffer(b, np.float32).astype(np.float64) for b in row]
    return arrays + [None] * (3 - len(arrays))

x_new, y_new, f_new = lawn_grass(REFERENCE_DB)
x_old, y_old, _ = lawn_grass(REPO / "resources" / "tetracorder_rules_references.db")
f_old = np.interp(x_old, x_new, f_new)      # Beckman FWHM onto the ASCII grid

for label, (x, y, f) in {"splib06a copy": (x_new, y_new, f_new),
                         "old ASCII copy": (x_old, y_old, f_old)}.items():
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(f)
    p = conv.convolve(x[keep], f[keep], y[keep])
    ok = valid_bands & np.isfinite(native) & np.isfinite(p)
    d = np.abs(p - native)
    worst = np.argmax(np.where(ok, d, 0))
    print(f"{label:15s} {keep.sum():4d} ch  max|diff| {d[ok].max():.2e} at "
          f"{wavelengths[worst]:.4f} um  (python NaN where native valid: "
          f"{int((valid_bands & np.isfinite(native) & np.isnan(p)).sum())})")
    
#%% hypothesis: native 7644 was convolved from a cubic-spline 3961-ch version
from scipy.interpolate import CubicSpline

con = sqlite3.connect(REFERENCE_DB)
gx, gf = (np.frombuffer(b, np.float32).astype(np.float64) for b in con.execute(
    "SELECT XData, FData FROM Samples JOIN Spectra USING (SampleID) "
    "WHERE Library = 'splib06b' AND ConvolvedRecord = 2568").fetchone())
con.close()
print(f"3961 grid from jarosite-K: {gx.size} ch, {gx[0]:.4f}-{gx[-1]:.4f} um")

keep = np.isfinite(x_new) & np.isfinite(y_new)
inside = (gx >= x_new[keep][0]) & (gx <= x_new[keep][-1])
for bc in ("not-a-knot", "natural"):
    y_csp = CubicSpline(x_new[keep], y_new[keep], bc_type=bc)(gx[inside])
    p = conv.convolve(gx[inside], gf[inside], y_csp)
    ok = valid_bands & np.isfinite(native) & np.isfinite(p)
    d = np.abs(p - native)
    worst = np.argmax(np.where(ok, d, 0))
    print(f"csp3961 ({bc:10s}) max|diff| {d[ok].max():.2e} at {wavelengths[worst]:.4f} um")
    
#%% benchmark _window_mean: old (span + mask) vs new (channel loop) on real features
import json, time
import numpy as np
from tetracorderp.tetracorder_ops import _wtochbin, prepare_rules

N_FEATURES = 20          # linear features to test (distinct window pairs)
REPEATS = 3              # timing repeats per call, best taken
BLOCK = slice(0, 100)    # cube lines to use

def window_mean_old(values, wav, mask):
    ok = mask & np.isfinite(values)
    n = ok.sum(axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        refl = (np.where(ok, values, np.float32(0)).sum(axis=-1, dtype=np.float32)
                / n.astype(np.float32))
        wl = (np.where(ok, wav, np.float32(0)).sum(axis=-1, dtype=np.float32)
              / n.astype(np.float32))
    return refl.astype(np.float32), wl.astype(np.float32), n

def window_mean_new(values, wav, mask):
    shape = values.shape[:-1]
    refl_sum = np.zeros(shape, dtype=np.float32)
    wl_sum = np.zeros(shape, dtype=np.float32)
    n = np.zeros(shape, dtype=np.intp)
    for c in np.flatnonzero(mask):
        v = values[..., c]
        good = np.isfinite(v)
        refl_sum += np.where(good, v, np.float32(0))
        wl_sum += np.where(good, np.float32(wav[c]), np.float32(0))
        n += good
    with np.errstate(invalid="ignore", divide="ignore"):
        count = n.astype(np.float32)
        return refl_sum / count, wl_sum / count, n

def best_time(fn, *args):
    t = []
    for _ in range(REPEATS):
        t0 = time.perf_counter(); out = fn(*args); t.append(time.perf_counter() - t0)
    return min(t), out

def ulps(a, b):
    """max float32 ulp distance where both finite"""
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return 0
    ia = a[ok].astype(np.float32).view(np.int32).astype(np.int64)
    ib = b[ok].astype(np.float32).view(np.int32).astype(np.int64)
    return int(np.abs(ia - ib).max())

# ---- inputs, prepared as linear_feature_continuum does --------------------
wl32 = np.asarray(wavelengths, dtype=np.float32)
valid = np.asarray(valid_bands, dtype=bool)
block = np.ascontiguousarray(np.asarray(cube[BLOCK], dtype=np.float32))

with open(RULES, encoding="utf-8") as fh:
    all_rules = prepare_rules(json.load(fh), mode=MODE)["rules"]

seen, features = set(), []
for mid, rule in all_rules.items():
    for f in rule["features"]:
        if f.get("role") == "not" or f["continuum"] != "linear":
            continue
        key = json.dumps(f["windows"])
        if key in seen:
            continue
        seen.add(key)
        features.append((mid, f["id"], f["windows"]))
    if len(features) >= N_FEATURES:
        break

# ---- run ------------------------------------------------------------------
print(f"block {block.shape}, {len(features)} features\n")
print(f"{'material':30s} {'feat':4s} {'n':>3s} {'kL':>3s} {'kR':>3s} "
      f"{'old ms':>8s} {'new ms':>8s} {'x':>5s}  n_eq  refl_ulp wl_ulp")
tot_old = tot_new = 0.0
for mid, fid, (lw, rw) in features:
    try:
        cl1, cl2 = _wtochbin(wl32, float(lw[0]), float(lw[1]))
        cr1, cr2 = _wtochbin(wl32, float(rw[0]), float(rw[1]))
    except ValueError:
        continue
    span = slice(cl1, cr2 + 1)
    wav = wl32[span]
    k = np.arange(cl1, cr2 + 1)
    left, right = k <= cl2, k >= cr1
    values = np.where(valid[span], block[..., span], np.nan).astype(np.float32)

    t_old = t_new = 0.0
    n_eq, r_ulp, w_ulp = True, 0, 0
    for mask in (left, right):
        to, (ro, wo, no) = best_time(window_mean_old, values, wav, mask)
        tn, (rn, wn, nn) = best_time(window_mean_new, values, wav, mask)
        t_old += to; t_new += tn
        n_eq &= np.array_equal(no, nn)
        r_ulp = max(r_ulp, ulps(ro, rn)); w_ulp = max(w_ulp, ulps(wo, wn))
    tot_old += t_old; tot_new += t_new
    print(f"{mid[:30]:30s} {fid:4s} {values.shape[-1]:3d} {left.sum():3d} {right.sum():3d} "
          f"{t_old*1e3:8.1f} {t_new*1e3:8.1f} {t_old/t_new:5.1f}  {str(n_eq):5s} "
          f"{r_ulp:8d} {w_ulp:6d}")

print(f"\ntotal: old {tot_old:.2f} s, new {tot_new:.2f} s, speed-up {tot_old/tot_new:.1f}x")
#%% benchmark characterise_feature sums: old (broadcast + products) vs new (matvec)
import time
import numpy as np
from tetracorderp.tetracorder_ops import linear_feature_continuum

N_FEATURES = 20
REPEATS = 3
BLOCK = slice(0, 100)
TINY = np.float32(0.1e-20)

def sums_old(rflibc, rfobsc):
    ok = np.isfinite(rfobsc) & np.isfinite(rflibc)
    n = ok.sum(axis=-1)
    drl = np.where(ok, rflibc, 0).astype(np.float64)
    dro = np.where(ok, rfobsc, 0).astype(np.float64)
    return (n, drl.sum(axis=-1), (drl * drl).sum(axis=-1), (dro * drl).sum(axis=-1),
            dro.sum(axis=-1), (dro * dro).sum(axis=-1))

def sums_new(rflibc, rfobsc):
    ref_ok = np.isfinite(rflibc)
    r = np.where(ref_ok, rflibc, 0).astype(np.float64)
    ok = np.isfinite(rfobsc) & ref_ok
    n = ok.sum(axis=-1)
    okf = ok.astype(np.float64)
    dro = np.where(ok, rfobsc, 0).astype(np.float64)
    return (n, okf @ r, okf @ (r * r), dro @ r,
            dro.sum(axis=-1), np.einsum("...i,...i->...", dro, dro))

def finish(sums, rflibc, ext):
    """downstream arithmetic of characterise_feature, unchanged"""
    n, suml, sumll, sumol, sumo, sumoo = sums
    dxn = np.where(n > 0, n, 1).astype(np.float64)
    top = (sumol - sumo * suml / dxn).astype(np.float32)
    bottom = (sumll - suml * suml / dxn).astype(np.float32)
    botm2 = (sumoo - sumo * sumo / dxn).astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        slope = np.where(np.abs(bottom) < TINY, 0, top / bottom).astype(np.float32)
        xk = ((1.0 - slope) / slope).astype(np.float32)
        xk1 = (xk + 1.0).astype(np.float32)
        depth = (1.0 - ((rflibc[ext] + xk) / xk1).astype(np.float32)).astype(np.float32)
        bprime = np.where(np.abs(botm2) < TINY, 0, top / botm2).astype(np.float32)
        fit = np.sqrt(np.abs(slope * bprime)).astype(np.float32)
    return dict(top=top, bottom=bottom, botm2=botm2, fit=fit, depth=depth)

def best_time(fn, *args):
    t = []
    for _ in range(REPEATS):
        t0 = time.perf_counter(); out = fn(*args); t.append(time.perf_counter() - t0)
    return min(t), out

def ulps(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return 0
    ia = a[ok].astype(np.float32).view(np.int32).astype(np.int64)
    ib = b[ok].astype(np.float32).view(np.int32).astype(np.int64)
    return int(np.abs(ia - ib).max())

def nan_mismatch(a, b):
    return int((np.isnan(a) != np.isnan(b)).sum())

# ---- inputs ----------------------------------------------------------------
ev = run.evaluator
wl32 = np.asarray(wavelengths, dtype=np.float32)
block = np.ascontiguousarray(np.asarray(cube[BLOCK], dtype=np.float32))

features, seen = [], set()
for mid, rule in ev.rules.items():
    if mid not in ev.rule_spectra:
        continue
    for f in rule["features"]:
        if f.get("role") == "not" or f["continuum"] != "linear":
            continue
        key = (mid, f["id"])
        if key in seen:
            continue
        seen.add(key)
        features.append((mid, f["id"], f["windows"]))
    if len(features) >= N_FEATURES:
        break

# ---- run -------------------------------------------------------------------
print(f"block {block.shape}, {len(features)} features\n")
print(f"{'material':28s} {'feat':4s} {'n':>3s} {'old ms':>8s} {'new ms':>8s} {'x':>5s}  "
      f"n_eq  sum_rel   top bot bm2 fit dep  nan")
tot_old = tot_new = 0.0
for mid, fid, (lw, rw) in features:
    reference = ev.rule_spectra[mid]
    valid = np.asarray(valid_bands, bool) & np.isfinite(np.asarray(reference, float))
    try:
        ref_c = linear_feature_continuum(reference, wl32, lw, rw, valid_bands=valid)
        tgt_c = linear_feature_continuum(block, wl32, lw, rw, valid_bands=valid)
    except ValueError:
        continue
    rflibc = np.asarray(np.ma.filled(np.ma.asarray(ref_c.continuum_removed, np.float32), np.nan))
    rfobsc = np.asarray(np.ma.filled(np.ma.asarray(tgt_c.continuum_removed, np.float32), np.nan))

    inner = np.flatnonzero(ref_c.interior & np.isfinite(rflibc))
    if inner.size == 0:
        continue
    minch = inner[np.argmin(rflibc[inner])]
    maxch = inner[np.argmax(rflibc[inner])]
    ext = maxch if (rflibc[maxch] - 1.0) > (1.0 - rflibc[minch]) else minch

    to, so = best_time(sums_old, rflibc, rfobsc)
    tn, sn = best_time(sums_new, rflibc, rfobsc)
    tot_old += to; tot_new += tn

    n_eq = np.array_equal(so[0], sn[0])
    rel = max(float(np.nanmax(np.abs(a - b) / np.maximum(np.abs(a), 1e-300)))
              for a, b in zip(so[1:], sn[1:]))
    fo, fn_ = finish(so, rflibc, ext), finish(sn, rflibc, ext)
    u = {k: ulps(fo[k], fn_[k]) for k in fo}
    nm = sum(nan_mismatch(fo[k], fn_[k]) for k in fo)

    print(f"{mid[:28]:28s} {fid:4s} {rflibc.size:3d} {to*1e3:8.1f} {tn*1e3:8.1f} "
          f"{to/tn:5.1f}  {str(n_eq):5s} {rel:8.1e} {u['top']:4d} {u['bottom']:3d} "
          f"{u['botm2']:3d} {u['fit']:3d} {u['depth']:3d} {nm:4d}")

print(f"\ntotal sums: old {tot_old:.2f} s, new {tot_new:.2f} s, speed-up {tot_old/tot_new:.1f}x")