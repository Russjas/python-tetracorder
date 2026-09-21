#%%
"""
Run the python-tetracorder package end to end on the native cuprite95 test
inputs and compare with the native Tetracorder 6.00 run.

This script only PREPARES inputs:
    cube          native int16 cube, (DN + offset) * scale, deleted DN -> NaN
    wavelengths   native wavelength record (s06av95a record 6)
    valid bands   native [DELETPTS]
    references    native SPECpr-convolved s06av95a / r06av95a records

and hands them to the package. Everything else is the package itself:
GroupEvaluator.evaluate() (groups, group 0, NOGROUP0, cases) and
MaterialEvaluator (features, tests, weights, constraints, NOTs).

The only substitution is where reference spectra come from: a
MaterialEvaluator subclass returns the native convolved spectra from
_prepare_rule_spectra / _prepare_rratio_spectra instead of reading the
SQLite DB and convolving. No package function is patched.
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
R06_AV95 = first_existing(
    WSL_ROOT / "sl1" / "usgs" / "rlib06" / "r06av95a",
    TET_ROOT / "sl1" / "usgs" / "rlib06" / "r06av95a",
    REPO / "scratch" / "r06av95a",
)
LIBRARY_FILES = {
    "[splib06]": S06_AV95, "[sprlb06]": R06_AV95,     # rule SMALL libraries
    "splib06b": S06_AV95, "sprlb06b": R06_AV95,       # config.RRATIO_REFERENCES
}
print("s06av95a:", S06_AV95)
print("r06av95a:", R06_AV95)


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


# =============================================================================
# THE ONLY SUBSTITUTION: reference spectra come from the native libraries
# =============================================================================

class NativeReferenceMaterialEvaluator(MaterialEvaluator):
    """MaterialEvaluator whose references are native SPECpr-convolved spectra."""

    def _prepare_rule_spectra(self, references):
        return native_rule_spectra(self.rules)

    def _prepare_rratio_spectra(self, references):
        return {name: read_specpr(LIBRARY_FILES[lib], rec)[1]
                for name, (lib, rec) in RRATIO_REFERENCES.items()}


class NativeReferenceGroupEvaluator(GroupEvaluator):
    """GroupEvaluator.__init__ with the evaluator above; evaluate() unchanged."""

    def __init__(self, target_spectra, target_wavelengths, mode="default",
                 target_valid_bands=None, rules_file=None,
                 disabled_groups=None, disabled_cases=None, disabled_materials=None,
                 temperature=None, pressure=None):
        self.target = target_spectra
        self.wavelengths = target_wavelengths
        self.fwhm = None
        self.mode = mode
        self.valid_bands = target_valid_bands
        if self.valid_bands is None:
            self.valid_bands = np.ones(self.wavelengths.shape, dtype=bool)
        self.reference_file = None
        self.rules_file = rules_file
        self.disabled_groups = {int(g) for g in (disabled_groups or ())}
        self.disabled_cases = {int(c) for c in (disabled_cases or ())}
        self.disabled_materials = set(disabled_materials or ())
        self.temperature = temperature
        self.pressure = pressure
        self.evaluator = NativeReferenceMaterialEvaluator(
            self.rules_file,
            {},                              # references come from the libraries
            self.wavelengths,
            self.fwhm,
            target_valid_bands=self.valid_bands,
            mode=mode,
            disabled_materials=self.disabled_materials,
        )
        self.disabled_materials |= self._resolve_physical_and_disable()
        self.evaluator.disabled_materials = self.disabled_materials
        self.group_winners, self.case_winners = self.evaluate()


#%% =============================================================================
# RUN THE PACKAGE
# =============================================================================
TEMPERATURE, PRESSURE = native_conditions()
FORCED_GROUPS, FORCED_CASES = native_force_disabled()
print(f"data conditions: {TEMPERATURE} K, {PRESSURE} bar")
print(f"force-disabled: groups {sorted(FORCED_GROUPS)}, cases {sorted(FORCED_CASES)}")

t0 = time.perf_counter()
run = NativeReferenceGroupEvaluator(
    cube, wavelengths, mode=MODE,
    target_valid_bands=valid_bands, rules_file=str(RULES),
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
from tetracorderp.config import VARIABLE_PRESETS, OUTPUT_DIRS


def write_plain_jpgs(self, output_dir):
           from matplotlib import colormaps
           from matplotlib.image import imsave
           output_dir = Path(output_dir)
           output_dir.mkdir(parents=True, exist_ok=True)
           rules = self.evaluator.rules
           nl, ns = self.target.shape[:2]
           palette = np.vstack([colormaps[c].colors for c in ("tab20", "tab20b", "tab20c")])
           palette = palette[np.ptp(palette, axis=1) > 0.1]          # same colours as write_jpgs
           for kind, results in (("group", self.group_winners), ("case", self.case_winners)):
               for number, result in sorted(results.items()):
                   if result is None:
                       continue
                   winner = np.asarray(result["winner"], dtype=object)
                   slots = [m for m, r in rules.items() if r["kind"] == kind
                            and int(r["number"]) in ({number, 0} if kind == "group" else {number})]
                   fd = np.ma.filled(np.ma.asarray(result["fit_depth"], dtype=np.float32), 0.0)
                   rgb = np.zeros((nl, ns, 3))
                   for mid in set(winner[winner != None]):  # noqa: E711
                       won = winner == mid
                       lines = [l.split("\\#")[0].strip() for l in rules[mid]["output_raw"].splitlines()]
                       lines = [l for l in lines if l]
                       dn, value = re.match(r"\d+\s+DN\s+(\d+)\s*=\s*(\S+)", lines[2]).groups()
                       value = VARIABLE_PRESETS[self.mode].get(value, value)
                       d = np.clip(np.floor(fd[won] * (int(dn) / float(value)) + 0.5), 0, 255)
                       w = (np.clip(d * 7.5 * 0.08 ** np.sqrt(d / 400.0), 0, 255) / 255.0)[:, None]
                       rgb[won] = palette[slots.index(mid) % len(palette)] * w      # gen.fd.gif.images gamma
                   imsave(output_dir / f"{kind}{number:02d}_{OUTPUT_DIRS[kind][number].split('.', 1)[1]}.jpg",
                          np.clip(rgb, 0, 1), pil_kwargs={"quality": 95})
write_plain_jpgs(run, "C:/Users/Hyperspectral/Documents/GitHub/python-tetracorder/scratch/write_jpgs_native_convolved_plain")

#%%
import re, numpy as np
from pathlib import Path
from matplotlib.image import imsave
from src.config import VARIABLE_PRESETS

support = Path(r"\\wsl.localhost\Ubuntu\home\hyperspectral\tetracorder-data\cuprite95\testrun1\cmds.color.support")
master = "\n".join(l for l in (support / "davinci.master.colors").read_text().splitlines()
                   if not l.lstrip().startswith("#"))
colours = {n: np.array([float(r), float(g), float(b)]) / 255.0 for r, g, b, n in re.findall(
    r"argb\[1,1,1\]\s*=\s*([\d.]+)\s*;\s*argb\[1,1,2\]\s*=\s*([\d.]+)\s*;\s*argb\[1,1,3\]\s*=\s*([\d.]+)\s*;\s*(c_\w+)\s*=\s*argb",
    master)}
clark, c = {}, None
for line in (support / "davinci.make.2micron-mins-detail2.a").read_text().splitlines():
    line = line.split("#")[0]
    if m := re.match(r"\s*c\s*=\s*(c_\w+)", line):
        c = colours[m.group(1)]
    if m := re.search(r'group\.2um/([^"]+?)\.(fd\.gif|depth\.gz)"', line):
        clark[m.group(1)] = (c, m.group(2))

rules = run.evaluator.rules
result = run.group_winners[2]
winner = np.asarray(result["winner"], dtype=object)
planes = {"fd.gif": np.ma.filled(np.ma.asarray(result["fit_depth"], dtype=np.float32), 0.0),
          "depth.gz": np.ma.filled(np.ma.asarray(result["depth"], dtype=np.float32), 0.0)}

rgb = np.zeros(winner.shape + (3,))
for mid in set(winner[winner != None]):  # noqa: E711
    lines = [l.split("\\#")[0].strip() for l in rules[mid]["output_raw"].splitlines()]
    lines = [l for l in lines if l]
    entry = clark.get(lines[1].split()[0])
    if entry is None:
        continue
    colour, source = entry
    won = winner == mid
    dn, value = re.match(r"\d+\s+DN\s+(\d+)\s*=\s*(\S+)", lines[2]).groups()
    value = VARIABLE_PRESETS[run.mode].get(value, value)
    d = np.clip(np.floor(planes[source][won] * (int(dn) / float(value)) + 0.5), 0, 255)
    if source == "fd.gif":
        d = np.clip(d * 7.5 * 0.08 ** np.sqrt(d / 400.0), 0, 255)      # gen.fd.gif.images gamma
    rgb[won] = colour * (d / 255.0)[:, None]                         # xcolor + a*(c/255)
imsave("group02_clark-detail2-emulatedbug.jpg", np.clip(rgb, 0, 1), pil_kwargs={"quality": 95})

#%%
import numpy as np
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.image import imread

run_dir = Path(r"\\wsl.localhost\Ubuntu\home\hyperspectral\tetracorder-data\cuprite95\testrun1")
nl, ns = rgb.shape[:2]
scene = imread(run_dir / "base-image" / "color-visRGB.jpg")
native = imread(run_dir / "color.results" /
                "cuprite95-aviris_t6.00geo2_color-results_2micron-minerals-detail2.png")[:, :ns, :3]

fig = Figure(figsize=(3 * 4 * ns / nl + 0.6, 4 + 0.5), facecolor="black")
axes = fig.subplots(1, 3, gridspec_kw={"wspace": 0.03})
for ax, image, title in zip(axes, (scene, native, np.clip(rgb, 0, 1)),
                            ("AVIRIS 1995, Cuprite NV", "Tetracorder 6.00 (ratfor)", "python-tetracorder")):
    ax.imshow(image, interpolation="nearest")
    ax.set_title(title, color="white", fontsize=11)
    ax.set_axis_off()
fig.text(0.5, 0.06, "2 µm minerals, group 2, coloured as Tetracorder davinci.make.2micron-mins-detail2.a",
         ha="center", color="0.7", fontsize=8)
fig.savefig("cuprite95_group2_comparison.jpg", dpi=600, bbox_inches="tight", facecolor="black",
            pil_kwargs={"quality": 95})