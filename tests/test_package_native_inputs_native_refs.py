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
# PACKAGE  (reloaded so edits to src are picked up in the same kernel)
# =============================================================================

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import src.config                                              # noqa: E402
import src.tetracorder_ops                                     # noqa: E402
import src.material_evaluation                                 # noqa: E402
import src.group_evaluator                                     # noqa: E402

for module in (src.config, src.tetracorder_ops,
               src.material_evaluation, src.group_evaluator):
    importlib.reload(module)

from src.config import RRATIO_REFERENCES                       # noqa: E402
from src.material_evaluation import MaterialEvaluator          # noqa: E402
from src.group_evaluator import GroupEvaluator                 # noqa: E402


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
import importlib, src.tetracorder_ops, src.material_evaluation
importlib.reload(src.tetracorder_ops)
importlib.reload(src.material_evaluation)
from src.config import NOGROUP0                                # noqa: E402


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
for label, rows in (("case 1", compare_set("case", 1, run.case_winners.get(1),
                     [m for m, r in rules.items() if r["kind"] == "case" and int(r["number"]) == 1])[0]),):
    for r in rows:
        print(f"{r['material']:22s} fit  n={r['fit']['n']:7d} exact={r['fit']['exact']:.3f} "
              f"<=1DN={r['fit']['within1']:.4f} median={r['fit']['median']:+.0f} max={r['fit']['max']}")
        print(f"{'':22s} dep  n={r['depth']['n']:7d} exact={r['depth']['exact']:.3f} "
              f"<=1DN={r['depth']['within1']:.4f} median={r['depth']['median']:+.0f} max={r['depth']['max']}")
    
#%%
"""
Literal port of nvres.r + bandmp.r (+ the nvres1mat material step), run against
the package on case-1 pixels, with native's own output DN alongside.

Run in the SAME kernel after test_package_native_inputs_All_groups.py (uses
cube, wavelengths, valid_bands, run, rules, read_specpr, LIBRARY_FILES,
RRATIO_REFERENCES, output_dirs, read_native_byte_gz, output_info, depth_scales).

Everything below follows the Ratfor line by line in float32: no masked arrays,
no vectorisation, sequential accumulation exactly as the do-loops run.
"""
import numpy as np

CASE = 1
N_WORST = 6          # pixels with the largest |python - native| fit DN
N_TYPICAL = 4        # pixels where the two agree, as a control
DELPT = -1.23e34

f4 = np.float32
ev = run.evaluator
resolved = run.case_winners[CASE]


def deleted(x):
    return (not np.isfinite(x)) or x <= DELPT / 2.0


def wtochbin(wav, w1, w2):
    """specpr wtochbin: first channel >= w1, last channel <= w2 (1-based -> 0-based)."""
    first = int(np.argmax(wav >= w1))
    last = int(len(wav) - 1 - np.argmax(wav[::-1] <= w2))
    return first, last


def bdmset(wav, ref, cl1, cl2, cr1, cr2):
    """Reference continuum removal + band min/max, as reflsetup does once at setup."""
    sl = sw = 0.0
    n = 0
    for i in range(cl1, cl2 + 1):
        if deleted(ref[i]):
            continue
        sl += float(ref[i]); sw += float(wav[i]); n += 1
    avlc, avwlc = f4(sl / n), f4(sw / n)

    sr = sw = 0.0
    n = 0
    for i in range(cr1, cr2 + 1):
        if deleted(ref[i]):
            continue
        sr += float(ref[i]); sw += float(wav[i]); n += 1
    avrc, avwrc = f4(sr / n), f4(sw / n)

    a = f4((avrc - avlc) / (avwrc - avwlc))
    b = f4(avrc - a * avwrc)

    rflibc = np.full(len(wav), DELPT, dtype=np.float32)
    for i in range(cl1, cr2 + 1):
        if deleted(ref[i]):
            continue
        contin = f4(a * wav[i] + b)
        if abs(contin) > 0.1e-20:
            rflibc[i] = f4(ref[i] / contin)

    minch = maxch = -1
    rmin, rmax = None, None
    for i in range(cl2 + 1, cr1):                  # interior only
        if deleted(rflibc[i]):
            continue
        if rmin is None or rflibc[i] < rmin:
            rmin, minch = rflibc[i], i
        if rmax is None or rflibc[i] > rmax:
            rmax, maxch = rflibc[i], i
    nftype = -1 if (rmax - 1.0) > (1.0 - rmin) else 1
    return rflibc, minch, maxch, nftype


def nvres(wav, rflibc, rfobs, vegspec, cl1, cl2, cr1, cr2, minch, maxch, nftype):
    """nvres.r then bandmp.r with imgflg = 0."""
    out = {}
    rfobs = rfobs.copy()
    vegspec = vegspec.copy()

    def window(c1, c2):
        su = sv = 0.0
        n = 0
        for i in range(c1, c2 + 1):
            if (deleted(wav[i]) or deleted(rflibc[i])
                    or deleted(vegspec[i]) or deleted(rfobs[i])):
                rfobs[i] = DELPT
                continue
            su += float(rfobs[i]); sv += float(vegspec[i]); n += 1
        if n < 1:
            return None, None
        return f4(su / n), f4(sv / n)

    avlcu, avlcv = window(cl1, cl2)
    avrcu, avrcv = window(cr1, cr2)
    if avlcu is None or avrcu is None:
        return None
    out.update(avlcu=avlcu, avlcv=avlcv, avrcu=avrcu, avrcv=avrcv)

    atop = f4(avlcu * avrcv - avrcu * avlcv)
    abot = f4(avrcu - avlcu)
    if abs(abot) < 0.1e-10:
        return None
    a = f4(atop / abot)
    b = f4((avlcv + a) / avlcu)
    if abs(b) < 0.1e-10:
        return None
    out.update(a=a, b=b)

    for i in range(cl1, cr2 + 1):                  # ratio, in place as the Ratfor does
        if deleted(vegspec[i]) or abs(vegspec[i]) < 0.1e-10:
            rfobs[i] = DELPT
            continue
        if deleted(rfobs[i]):
            vegspec[i] = DELPT
            continue
        vegspec[i] = f4((vegspec[i] + a) / b)
        rfobs[i] = f4(rfobs[i] / vegspec[i])

    # ---- bandmp on the ratioed spectrum -------------------------------------
    def contwin(c1, c2):
        s = sw = 0.0
        n = 0
        for i in range(c1, c2 + 1):
            if deleted(wav[i]) or deleted(rflibc[i]) or deleted(rfobs[i]):
                continue
            s += float(rfobs[i]); sw += float(wav[i]); n += 1
        if n < 1:
            return None, None
        return f4(s / n), f4(sw / n)

    avlc, avwlc = contwin(cl1, cl2)
    avrc, avwrc = contwin(cr1, cr2)
    if avlc is None or avrc is None or avlc < 0.1e-20 or avrc < 0.1e-20:
        return None
    aa = f4((avrc - avlc) / (avwrc - avwlc))
    bb = f4(avrc - aa * avwrc)

    rfobsc = np.full(len(wav), DELPT, dtype=np.float32)
    for i in range(cl1, cr2 + 1):                  # imgflg = 0 -> il, ir = cl1, cr2
        if deleted(wav[i]) or deleted(rflibc[i]) or deleted(rfobs[i]):
            continue
        contin = f4(aa * wav[i] + bb)
        if abs(contin) > 0.1e-20:
            rfobsc[i] = f4(rfobs[i] / contin)

    ext = minch if nftype == 1 else maxch
    conref = f4(aa * wav[ext] + bb)
    out["conref"] = conref

    suml = sumll = sumol = sumo = sumoo = 0.0
    n = 0
    for i in range(cl1, cr2 + 1):
        if deleted(rfobsc[i]):
            continue
        n += 1
        drl, dro = float(rflibc[i]), float(rfobsc[i])
        suml += drl; sumll += drl * drl; sumol += dro * drl
        sumo += dro; sumoo += dro * dro
    if n == 0:
        return None
    out["n"] = n

    top = f4(sumol - sumo * suml / n)
    bottom = f4(sumll - suml * suml / n)
    bslope = f4(0.0) if abs(bottom) < 0.1e-20 else f4(top / bottom)
    if abs(bslope) < 0.1e-20:
        return None
    xk = f4((1.0 - bslope) / bslope)
    xk1 = f4(xk + 1.0)
    if abs(xk1) < 0.1e-20:
        return None
    bd = f4(1.0 - f4((rflibc[ext] + xk) / xk1))
    botm2 = f4(sumoo - sumo * sumo / n)
    bprime = f4(0.0) if abs(botm2) < 0.1e-20 else f4(top / botm2)
    fit = f4(np.sqrt(abs(bslope * bprime)))
    out.update(xk=xk, bd=bd, fit=fit)

    xndvivegspec = f4((avrcv - avlcv) / (avrcv + avlcv))
    xndviobs = f4((avrcu - avlcu) / (avrcu + avlcu))
    if abs(xndviobs) <= 0.1e-10:
        return None
    xfactor = f4(xndvivegspec / xndviobs)
    out.update(xfactor=xfactor, bdnorm=f4(xfactor * bd), nftype=nftype)
    return out


def native_material(rule, r):
    """nvres1mat: ct, presence, weighted sums (one feature, dln = 1), FITALL."""
    if r is None:
        return 0.0, 0.0
    tests = {t["test"]: t["values"] for t in rule["features"][0].get("tests", [])}
    ct = tests.get("ct")
    ct = ct[0] if isinstance(ct, (list, tuple)) else ct
    fit, bdnorm = float(r["fit"]), float(r["bdnorm"])
    if ct is not None and r["conref"] < ct:
        return 0.0, 0.0
    xfeat = float(r["nftype"])
    xx = 0.0 if bdnorm / xfeat <= 0.1e-5 else 1.0
    ofit, odepth = fit * xx, bdnorm * xfeat
    for c in rule.get("constraints", []):
        if c["test"] != "FITALL>":
            continue
        z1, z2 = c["values"]
        if ofit < z1:
            return 0.0, 0.0
        if ofit < z2:
            f = (ofit - z1) / (z2 - z1)
            ofit, odepth = ofit * f, odepth * f
    return ofit, odepth


# ---------------------------------------------------------------- pixels
members = [m for m, r in rules.items() if r["kind"] == "case" and int(r["number"]) == CASE]
winner = resolved["winner"]
picks = []
for mid in members:
    basename, _ = output_info(rules[mid])
    path = output_dirs[("case", CASE)] / f"{basename}.fit.gz"
    if not path.exists():
        continue
    n_fit = read_native_byte_gz(path)
    p_fit = nint_uint8(np.where(winner == mid, resolved["fit"], 0.0), 255.0)
    both = (n_fit > 0) & (p_fit > 0)
    diff = np.where(both, np.abs(p_fit.astype(int) - n_fit.astype(int)), -1)
    order = np.argsort(diff, axis=None)[::-1]
    worst = [np.unravel_index(k, diff.shape) for k in order[:N_WORST]]
    same = [np.unravel_index(k, diff.shape)
            for k in np.flatnonzero((diff == 0).ravel())[:N_TYPICAL]]
    picks.append((mid, basename, n_fit, worst + same))

for mid, basename, n_fit, pixels in picks:
    rule = rules[mid]
    feature = rule["features"][0]
    W = feature["windows"]
    ref = ev.rule_spectra[mid]
    ratio_name = rule["preratio"].split()[1]
    veg = ev.rratio_spectra[ratio_name]
    depth_scale = depth_scales.get(("case", CASE, basename), 510.0)

    n_depth = read_native_byte_gz(output_dirs[("case", CASE)] / f"{basename}.depth.gz")
    cl1, cl2 = wtochbin(wavelengths, W[0][0], W[0][1])
    cr1, cr2 = wtochbin(wavelengths, W[1][0], W[1][1])
    rflibc, minch, maxch, nftype = bdmset(wavelengths, ref, cl1, cl2, cr1, cr2)

    print("\n" + "=" * 78)
    print(f"{mid}  ({basename})  ratio {ratio_name}")
    print(f"  channels {cl1+1}-{cl2+1} .. {cr1+1}-{cr2+1}  band ext "
          f"{(minch if nftype == 1 else maxch)+1}  nftype {nftype}")
    print("=" * 78)

    for row, col in pixels:
        obs = np.asarray(cube[row, col], dtype=np.float32)
        obs = np.where(np.isfinite(obs) & valid_bands, obs, DELPT).astype(np.float32)
        vegpx = np.where(np.isfinite(veg) & valid_bands, veg, DELPT).astype(np.float32)

        r = nvres(wavelengths.astype(np.float32), rflibc, obs, vegpx,
                  cl1, cl2, cr1, cr2, minch, maxch, nftype)
        ofit, odepth = native_material(rule, r)

        p_fit = float(resolved["fit"][row, col]) if winner[row, col] == mid else 0.0
        p_depth = float(resolved["depth"][row, col]) if winner[row, col] == mid else 0.0

        print(f"\n  line {row} sample {col}   native DN fit {n_fit[row, col]:3d} "
              f"depth {n_depth[row, col]:3d}")
        if r is None:
            print("    port: deleted")
        else:
            print(f"    port  avlcu {r['avlcu']:.6f} avlcv {r['avlcv']:.6f} "
                  f"avrcu {r['avrcu']:.6f} avrcv {r['avrcv']:.6f}")
            print(f"    port  a {r['a']:.6f} b {r['b']:.6f} conref {r['conref']:.6f} "
                  f"n {r['n']} bd {r['bd']:.6f} xfactor {r['xfactor']:.6f}")
            print(f"    port  fit {r['fit']:.6f} bdnorm {r['bdnorm']:.6f}  "
                  f"-> material fit {ofit:.6f} depth {odepth:.6f}  "
                  f"DN fit {int(np.floor(ofit*255+0.5))} "
                  f"depth {int(np.floor(abs(odepth)*depth_scale+0.5))}")
        print(f"    pkg   fit {p_fit:.6f} depth {p_depth:.6f}  "
              f"DN fit {int(np.floor(p_fit*255+0.5))} "
              f"depth {int(np.floor(abs(p_depth)*depth_scale+0.5))}")
    
#%%
import inspect, src.material_evaluation as me
print(me.__file__)
s = inspect.getsource(me.resolve_feature_tests)
i = s.index('if "r*bd>" in tests_by_name')
print(s[i:i+700])
#%%
g1 = run.group_winners[1]
mid = "fe3+bearing1"
basename, _ = output_info(rules[mid])
n_fit = read_native_byte_gz(output_dirs[("group", 1)] / f"{basename}.fit.gz")
mine = g1["winner"] == mid
p_fit = nint_uint8(np.where(mine, g1["fit"], 0.0), 255.0)

only = (p_fit > 0) & (n_fit == 0)
both = (p_fit > 0) & (n_fit > 0)
print("python-only", int(only.sum()), " both", int(both.sum()))
for label, m in (("python-only", only), ("both", both)):
    raw = g1["fit"][m]
    print(f"{label:12s} raw fit  min {raw.min():.4f}  median {np.median(raw):.4f}  "
          f"max {raw.max():.4f}   DN median {int(np.median(p_fit[m]))}")
#%%
import src.material_evaluation as me
mid = "fe3+bearing1"
pix = np.argwhere(only)[:4]
spectra = cube[pix[:, 0], pix[:, 1]][:, None, :]
ev = run.evaluator
rule = ev.rules[mid]
feats = {f["id"]: f for f in rule["features"]}
res = ev._resolve_positive_features(feats, ev.rule_spectra[mid], spectra)
w = me.resolve_feature_weights(feats, res)
flat = lambda a: np.ma.filled(np.ma.asarray(a, float), np.nan).ravel()

for fid, r in res.items():
    print(f"{fid} {feats[fid]['role']:10s} {r['status']:8s} {str(r['polarity']):10s} "
          f"w={w[fid]:.3f}")
    print("     raw fit", flat(r["raw_fit"]), " mod fit", flat(r["mod_fit"]))
    print("     raw dep", flat(r["raw_depth"]), " mod dep", flat(r["mod_depth"]))

mat = me.resolve_positive_material(feats, res)
con = me.resolve_material_constraints(mat, rule.get("constraints", []))
print("sum fit", flat(mat["fit"]), "-> after constraints", flat(con["fit"]))
#%%
members = [m for m, r in rules.items() if r["kind"] == "group" and int(r["number"]) == 1]
members += [m for m, r in rules.items() if r["kind"] == "group" and int(r["number"]) == 0]
any_native = np.zeros(cube.shape[:2], bool)
for m in members:
    bn, bits = output_info(rules[m])
    p = output_dirs[("group", 1)] / f"{bn}.fit.gz"
    if bn and bits == 8 and p.exists():
        any_native |= read_native_byte_gz(p) > 0

empty = only & ~any_native          # native found nothing at all in group 1
swap  = only & any_native           # native picked someone else
print("native empty", int(empty.sum()), "  swap", int(swap.sum()))
for label, m in (("empty", empty), ("swap", swap)):
    f = g1["fit"][m]
    print(f"{label:6s} final fit median {np.median(f):.4f} max {f.max():.4f}; "
          f"DN median {int(np.median(p_fit[m]))} max {int(p_fit[m].max())}")
ids, counts = np.unique(g1['winner'][swap].astype(str), return_counts=True)
#%%
vals = g1["fit"][empty]
idx = np.argwhere(empty)
order = np.argsort(vals)
pix = idx[order[[0, len(order)//4, len(order)//2, len(order)*3//4, -1]]]
spectra = cube[pix[:, 0], pix[:, 1]][:, None, :]

ff, st = me.fit_feature(ev.rule_spectra[mid], spectra, ev.wavelengths,
                        feats["f1a"]["windows"], continuum="linear",
                        valid_bands=ev.valid_bands)
print("status", st)
lc, rc = flat(ff["left_continuum"]), flat(ff["right_continuum"])
print("conref ", np.round(flat(ff["continuum"]), 4))
print("lc     ", np.round(lc, 4))
print("rc     ", np.round(rc, 4))
print("rc/lc  ", np.round(rc / lc, 4), "   (test: >0.9, full at 1.1)")
print("fit    ", np.round(flat(ff["fit"]), 4))
print("depth  ", np.round(flat(ff["depth"]), 5))
print("ext ch ", ff["extremum_index"] + 1, " polarity", ff["polarity"])
#%%
rows, cols = np.nonzero(empty)
print("lines", rows.min(), rows.max(), " samples", cols.min(), cols.max())
print("by line decile", np.histogram(rows, bins=10)[0])
print("by sample decile", np.histogram(cols, bins=10)[0])
    #%%
    
edge = np.zeros_like(empty)
edge[:, 560:] = True
p = np.argwhere(empty & edge)[0]
sp = cube[p[0], p[1]]
print("pixel", p, " nan", int(np.isnan(sp).sum()), "/", sp.size,
      " min", np.nanmin(sp), " max", np.nanmax(sp))
print(np.round(sp[:30], 4))
for mid2 in ("natrolite", "analcime"):
    w = run.group_winners[2]["winner"] == mid2
    r, c = np.nonzero(w)
    print(f"{mid2:10s} n={w.sum():4d}  samples {c.min()}-{c.max()}  "
          f"deciles {np.histogram(c, bins=10, range=(0, 614))[0]}")
    #%%
anyg = np.zeros(cube.shape[:2], bool)
for gnum in (2, 3, 4):
    d = output_dirs[("group", gnum)]
    for m, r in rules.items():
        if r["kind"] != "group":
            continue
        if int(r["number"]) not in (gnum, 0):
            continue
        bn, bits = output_info(r)
        p = d / f"{bn}.fit.gz"
        if bn and bits == 8 and p.exists():
            anyg |= read_native_byte_gz(p) > 0

print("native has some detection anywhere in groups 2-4:", int(anyg.sum()))
print("of the 5821 'empty' pixels, native silent everywhere too:",
      int((empty & ~anyg).sum()))
print("natrolite pixels where native is silent everywhere:",
      int(((run.group_winners[2]['winner'] == 'natrolite') & ~anyg).sum()))
#%%
d = output_dirs[("group", 1)]
have = {p.name[:-len(".fit.gz")] for p in d.glob("*.fit.gz")}
want = {}
for m, r in rules.items():
    if r["kind"] == "group" and int(r["number"]) in (1, 0):
        bn, bits = output_info(r)
        want[bn] = (m, bits)
print("files with no rule :", sorted(have - set(want)))
print("rules with no file :", sorted(set(want) - have))
print("of those, disabled :", sum(1 for bn in set(want) - have
                                  if want[bn][0] in run.disabled_materials))
win = np.full(cube.shape[:2], "", object)
for bn, (m, bits) in want.items():
    p = d / f"{bn}.fit.gz"
    if bits == 8 and p.exists():
        dn = read_native_byte_gz(p)
        win[dn > 0] = bn
print("native group-1 winner on the 5821:", 
      dict(zip(*np.unique(win[empty], return_counts=True))))
#%%
sel = []
for m, msk in (("fe3", empty), ("nat", run.group_winners[2]["winner"] == "natrolite"),
               ("ana", run.group_winners[2]["winner"] == "analcime")):
    idx = np.argwhere(msk)
    vals = (g1["fit"] if m == "fe3" else run.group_winners[2]["fit"])[msk]
    order = np.argsort(vals)
    take = idx[order[[len(order)//2, -1]]] if len(order) > 1 else idx
    sel += [(int(r), int(c), m) for r, c in take]
print(sel)
#%%
"""
Write the disputed pixels into a SPECpr file and generate the tetracorder
single-spectrum command script for them.

Run in the SAME kernel, after test_package_native_inputs_All_groups.py and
after the cell that defines `empty` (native silent in group 1 where Python
picks fe3+bearing1). Uses: run, cube, wavelengths, valid_bands, NATIVE_RUN,
LIBRARY_FILES, _physical_record, read_specpr, _LIBRARY_BYTES.

Native single-spectrum mode at diagnostic level 3 prints, per material and
feature, zfit / zdepth / zfd / dln / zcompf and a line per NOT feature - the
same quantities the Python side prints, for the same pixels.
"""
import struct
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------- settings
OUT_DIR = WSL_ROOT / "home" / "hyperspectral" 
SP_NAME = "pixels.sp"                     # SPECpr file of pixel spectra
SCRIPT_NAME = "cmds.single.pixels"        # tetracorder stdin script
FILE_LETTER = "v"                         # restart-file slot to assign it to
START_FILE = "cmds.start.t6.00a.single"   # check the name in the run directory

TEMPLATE_RECORD = 666                     # any 224-channel record: header donor
HEADER_RECORDS = 6                        # label + text records copied verbatim
DIAG_LEVEL = 3

DELETED = -1.23e34
RECORD_BYTES = 1536
TITLE_OFF, NCH_OFF, IRECNO_OFF, ITPNTR_OFF, DATA_OFF = 4, 80, 108, 112, 512
MAX_FIRST_RECORD_CHANNELS = 256

library = LIBRARY_FILES["[splib06]"]

# ---------------------------------------------------------------- pixels
g1 = run.group_winners[1]
g2 = run.group_winners[2]


def median_and_worst(mask, fit, label):
    """The median-fit pixel and the strongest one, as (line, sample, label)."""
    idx = np.argwhere(mask)
    if idx.size == 0:
        print(f"  no pixels for {label}")
        return []
    order = np.argsort(np.asarray(fit)[mask])
    take = idx[order[[len(order) // 2, -1]]] if len(order) > 1 else idx
    return [(int(r), int(c), label) for r, c in take]


PIXELS = (median_and_worst(empty, g1["fit"], "fe3+bearing1")
          + median_and_worst(g2["winner"] == "natrolite", g2["fit"], "natrolite")
          + median_and_worst(g2["winner"] == "analcime", g2["fit"], "analcime"))

if not PIXELS:
    raise SystemExit("no pixels selected")

for row, col, label in PIXELS:
    winner = g1["winner"][row, col] if label == "fe3+bearing1" else g2["winner"][row, col]
    fit = g1["fit"][row, col] if label == "fe3+bearing1" else g2["fit"][row, col]
    print(f"  line {row:4d} sample {col:4d}  {label:14s} python fit {fit:.4f} "
          f"(DN {int(np.floor(fit * 255 + 0.5))})  winner {winner}")

nbands = cube.shape[-1]
if nbands > MAX_FIRST_RECORD_CHANNELS:
    raise SystemExit(f"{nbands} channels needs continuation records; not implemented")

# ---------------------------------------------------------------- SPECpr file
template = _physical_record(library, TEMPLATE_RECORD)
out = bytearray()
for n in range(HEADER_RECORDS):                    # PDS label + text records
    out += _physical_record(library, n)

records = []
for i, (row, col, label) in enumerate(PIXELS):
    record_number = HEADER_RECORDS + i
    spectrum = np.asarray(cube[row, col], dtype=np.float64)
    values = np.where(np.isfinite(spectrum), spectrum, DELETED).astype(">f4")

    rec = bytearray(template)
    title = f"L{row} S{col} {label}"[:40].ljust(40)
    rec[TITLE_OFF:TITLE_OFF + 40] = title.encode("latin-1")
    struct.pack_into(">i", rec, NCH_OFF, nbands)
    struct.pack_into(">i", rec, IRECNO_OFF, record_number)
    struct.pack_into(">i", rec, ITPNTR_OFF, 0)
    rec[DATA_OFF:DATA_OFF + 4 * nbands] = values.tobytes()
    out += rec
    records.append((record_number, row, col, label))

sp_path = Path(OUT_DIR) / SP_NAME
sp_path.write_bytes(bytes(out))
print(f"\nwrote {sp_path}  ({len(PIXELS)} spectra, records "
      f"{records[0][0]}-{records[-1][0]}, {nbands} channels)")

# read back through the same reader the comparison uses
_LIBRARY_BYTES.pop(Path(sp_path), None)
for record_number, row, col, label in records:
    title, data = read_specpr(sp_path, record_number)
    same = np.array_equal(np.nan_to_num(data), np.nan_to_num(cube[row, col]))
    print(f"  rec {record_number}: {title!r:34s} finite={int(np.isfinite(data).sum()):3d} "
          f"round-trip={'ok' if same else 'MISMATCH'}")

# ---------------------------------------------------------------- run script
lines = [f"<{START_FILE}", "s"]
for record_number, row, col, label in records:
    lines += [
        f"{FILE_LETTER} {record_number}",       # file letter and record
        "",                                     # no min/max data thresholds
        str(DIAG_LEVEL),                        # full diagnostic output
        f"L{row} S{col} {label}",               # output comment, marks it in results
    ]
lines += ["x", "x", ""]

script_path = Path(OUT_DIR) / SCRIPT_NAME
with open(script_path, "w", newline="\n") as f:
    f.write("\n".join(lines))
print(f"wrote {script_path}")

# ---------------------------------------------------------------- restart patch
setup = f"""#!/bin/bash
# Point the restart file's {FILE_LETTER} slot at {SP_NAME}, in this directory.
set -e
here=$(cd "$(dirname "$0")" && pwd)
python3 - "$here" <<'PY'
import pathlib, re, sys
here = pathlib.Path(sys.argv[1])
path = here / "r1"
text = path.read_text()
text = re.sub(r"^i{FILE_LETTER}fl=.*$",
              "i{FILE_LETTER}fl=" + str(here / "{SP_NAME}").ljust(80),
              text, flags=re.M)
text = re.sub(r"^isavt=.*$",
              "isavt=      " + "pixels".rjust(8) + "  # file device letter {FILE_LETTER}",
              text, flags=re.M)
path.write_text(text)
print("patched", path)
PY
"""
setup_path = Path(OUT_DIR) / "setup_single.sh"
with open(setup_path, "w", newline="\n") as f:
    f.write(setup)
print(f"wrote {setup_path}")

print(f"""
In the container, on a COPY of the run directory - single-spectrum mode
rewrites history and results:

  cp -a {Path(OUT_DIR).name} singlerun && cd singlerun
  bash setup_single.sh
  /usr/local/bin/tetracorder6.00single r1 < {SCRIPT_NAME} 2>&1 | tee single.out

Then in `results`, each pixel appears under its comment (L<line> S<sample>),
followed by "FITS, DEPTHS, F*D before best fit selection" - native's own
zfit / zdepth / zfd / dln / zcompf per material and feature. The ones to read:

  material 105  fe3+bearing1      (feature 1 is the diagnostic 0.44-0.6435)
  material 326  natrolite
  material 316  analcime

If specpr complains that the record count of {SP_NAME} disagrees with the
restart file, answer c.
""")
#%%
idx = np.argwhere(empty)
row, col = idx[len(idx) // 2]
hdr = read_envi_header(str(cube_path) + ".hdr")
print(hdr, "deleted_dn", deleted_dn, "offset", dn_offset, "scale", dn_scale)

dtype = np.dtype({1: "u1", 2: "i2", 3: "i4", 4: "f4", 12: "u2"}[hdr["data_type"]])
dtype = dtype.newbyteorder(">" if hdr["byte_order"] == 1 else "<")
shape = {"bil": (hdr["lines"], hdr["bands"], hdr["samples"]),
         "bip": (hdr["lines"], hdr["samples"], hdr["bands"]),
         "bsq": (hdr["bands"], hdr["lines"], hdr["samples"])}[hdr["interleave"]]
raw = np.memmap(cube_path, dtype=dtype, mode="r", offset=hdr["offset"], shape=shape)

if hdr["interleave"] == "bil":
    dn = raw[row, :, col]
elif hdr["interleave"] == "bip":
    dn = raw[row, col, :]
else:
    dn = raw[:, row, col]
dn = np.asarray(dn)

print("pixel line", int(row), "sample", int(col))
print("raw DN  ", dn[:12])
print("== deleted:", int((dn == deleted_dn).sum()), "channels")
print("cube[]  ", np.round(cube[row, col][:12], 4))
print("expected", np.round((dn[:12].astype(np.float32) + dn_offset) * dn_scale, 4))
print("nan in cube:", int(np.isnan(cube[row, col]).sum()))
#%%
np.set_printoptions(precision=6, suppress=True, threshold=10000, linewidth=120)
print(repr(np.asarray(cube[626, 7])))
for gnum in (1, 2, 3, 4):
    d = output_dirs[("group", gnum)]
    hits = []
    for m, r in rules.items():
        if r["kind"] != "group" or int(r["number"]) not in (gnum, 0):
            continue
        bn, bits = output_info(r)
        p = d / f"{bn}.fit.gz" if bn else None
        if p is not None and bits == 8 and p.exists():
            v = int(read_native_byte_gz(p)[626, 7])
            if v:
                hits.append((bn, v))
    print(f"group {gnum}:", hits)
#%%
"""
Is native's group-1 winner a GROUP 0 material at the disputed pixels?

tp1all.r puts every group-0 material into the contest for each group not in
`nogroup0: 3 10 11 12`, but the "now that best is found, others are zero" loop
runs `do igroup = 1, nzgroup` - it never touches matgrp(*,0). So when a group-0
material wins group 1, all of group 1 is zeroed and group 1's raster is blank,
while groups 2 and 3 are unaffected. That is exactly what native shows on the
5821 pixels where Python maps fe3+bearing1.

Run in the SAME kernel as test_package_native_inputs_All_groups.py, after the
cell that defines `empty`. Uses: run, rules, cube, empty, output_dirs,
output_info, read_native_byte_gz, nint_uint8.
"""
import numpy as np

SAMPLE = 400          # disputed pixels to re-evaluate on the Python side
TOP = 12              # how many materials to list

g1 = run.group_winners[1]
group0_ids = [m for m, r in rules.items()
              if r["kind"] == "group" and int(r["number"]) == 0]
group1_ids = [m for m, r in rules.items()
              if r["kind"] == "group" and int(r["number"]) == 1]
print(f"group 0: {len(group0_ids)} materials, group 1: {len(group1_ids)}")

rows, cols = np.nonzero(empty)
print(f"disputed (native group-1 blank, python = fe3+bearing1): {rows.size}")

# ------------------------------------------------------------------ native
# What does native map at these pixels in the group-0 output directory?
g0_dir = output_dirs.get(("group", 0))
print(f"\nnative group-0 output dir: {g0_dir}")

native_g0 = {}
for mid in group0_ids:
    basename, bits = output_info(rules[mid])
    if g0_dir is None or basename is None or bits != 8:
        continue
    path = g0_dir / f"{basename}.fit.gz"
    if not path.exists():
        continue
    dn = read_native_byte_gz(path)
    hits = dn[rows, cols]
    if np.any(hits > 0):
        native_g0[mid] = (int((hits > 0).sum()), int(hits.max()),
                          float(np.median(hits[hits > 0])))

print(f"\n  native group-0 materials firing on the disputed pixels "
      f"({len(native_g0)} of {len(group0_ids)}):")
for mid, (n, mx, med) in sorted(native_g0.items(), key=lambda kv: -kv[1][0])[:TOP]:
    print(f"    {mid:44s} n={n:6d}  max DN {mx:3d}  median DN {med:5.1f}")
if not native_g0:
    print("    NONE - group 0 is not what is silencing group 1")

# ------------------------------------------------------------------ python
# Re-run every group-0 material on a sample of the disputed spectra.
take = np.linspace(0, rows.size - 1, min(SAMPLE, rows.size)).astype(int)
r, c = rows[take], cols[take]
spectra = np.asarray(cube[r, c])                     # (N, bands)

print(f"\n  python fits on {r.size} sampled disputed pixels:")
py_fit = {}
for mid in group0_ids + group1_ids:
    result = run.evaluator.evaluate(mid, spectra)
    if result is None:
        py_fit[mid] = None
        continue
    fit = np.ma.filled(np.ma.asarray(result["fit"]), 0.0).astype(float)
    py_fit[mid] = fit

dropped = [m for m, f in py_fit.items() if f is None]
print(f"    evaluate() returned None for {len(dropped)} materials"
      + (f": {dropped[:8]}{' ...' if len(dropped) > 8 else ''}" if dropped else ""))

live = {m: f for m, f in py_fit.items() if f is not None}
best_g0 = {m: f for m, f in live.items() if m in group0_ids}
fe3 = live.get("fe3+bearing1")

print(f"\n    fe3+bearing1 python fit on the sample: "
      f"median {np.median(fe3):.4f}  max {fe3.max():.4f}")
print("\n    group-0 materials by how often they BEAT fe3+bearing1 here:")
beat = {m: int((f > fe3).sum()) for m, f in best_g0.items()}
for mid, n in sorted(beat.items(), key=lambda kv: -kv[1])[:TOP]:
    f = best_g0[mid]
    print(f"    {mid:44s} beats on {n:4d}/{r.size}  median fit {np.median(f):.4f}"
          f"  max {f.max():.4f}")

any_beat = np.zeros(r.size, dtype=bool)
for f in best_g0.values():
    any_beat |= f > fe3
print(f"\n    ANY group-0 material beats fe3+bearing1 on "
      f"{int(any_beat.sum())}/{r.size} sampled pixels")

# ------------------------------------------------------------------ verdict
print("""
Reading it:
  native group-0 firing + python group-0 beating   -> ordering/None bug in
        resolve_group: the group-0 winner should have zeroed all of group 1.
  native group-0 firing + python group-0 silent    -> the group-0 material
        itself is mis-evaluated in Python; chase THAT material, not fe3+bearing1.
  native group-0 silent                            -> group 0 is not the cause;
        native is zeroing fe3+bearing1 inside tp1mat and we go back to bandmp.
""")
#%%
#%%
"""
Group 0, take 2: find the native rasters properly, then ask whether Python's
group 0 is alive at all.

cmds.start.t6.00a has no ==[DIRg0], so group-0 materials inherit ==[DIR]./ and
their .fit.gz sits in the run directory itself. native_output_dirs() only ever
makes ("group", N>=1) keys, so compare_set() has been looking for group-0
basenames inside the group-N directory and dropping all 48 as "without native
output". Locate them by basename instead.

Run in the SAME kernel, after the cell that defines `empty`. Uses: run, rules,
cube, empty, NATIVE_RUN, output_info, read_native_byte_gz.
"""
import numpy as np

SAMPLE = 2000          # random scene pixels for the Python side
TOP = 15

group0_ids = [m for m, r in rules.items()
              if r["kind"] == "group" and int(r["number"]) == 0]

# ------------------------------------------------------- locate the rasters
index = {}
for path in NATIVE_RUN.rglob("*.fit.gz"):
    index.setdefault(path.name[:-len(".fit.gz")], path)
print(f"indexed {len(index)} native .fit.gz files under {NATIVE_RUN}")

found, absent = {}, []
for mid in group0_ids:
    basename, bits = output_info(rules[mid])
    if basename is None or bits != 8:
        absent.append((mid, "no 8-bit output block"))
        continue
    path = index.get(basename)
    if path is None:
        absent.append((mid, f"no file for {basename!r}"))
        continue
    found[mid] = path

print(f"group 0: {len(found)}/{len(group0_ids)} have a native raster")
for mid, why in absent[:8]:
    print(f"    missing {mid:42s} ({why})")
if len(absent) > 8:
    print(f"    ... and {len(absent) - 8} more")
if found:
    print(f"    they live in: "
          f"{sorted({p.parent.relative_to(NATIVE_RUN).as_posix() or '.' for p in found.values()})}")

# ------------------------------------------------------- native, scene-wide
rows, cols = np.nonzero(empty)
print(f"\ndisputed pixels (native group-1 blank, python = fe3+bearing1): {rows.size}")
print(f"\n  {'material':44s} {'scene px':>9s} {'on disputed':>12s} {'max DN':>7s}")
native_counts = {}
for mid, path in found.items():
    dn = read_native_byte_gz(path)
    hits = dn[rows, cols]
    native_counts[mid] = (int((dn > 0).sum()), int((hits > 0).sum()), int(dn.max()))
for mid, (scene, disp, mx) in sorted(native_counts.items(), key=lambda kv: -kv[1][0])[:TOP]:
    print(f"  {mid:44s} {scene:9d} {disp:12d} {mx:7d}")

live_native = {m: v for m, v in native_counts.items() if v[0] > 0}
print(f"\n  native fires on {len(live_native)}/{len(found)} group-0 materials somewhere"
      f" in the scene; on the disputed pixels: "
      f"{sum(1 for v in native_counts.values() if v[1] > 0)}")

# ------------------------------------------------------- python, scene-wide
rng = np.random.default_rng(0)
lin = rng.choice(cube.shape[0] * cube.shape[1], size=SAMPLE, replace=False)
r, c = np.unravel_index(lin, cube.shape[:2])
spectra = np.asarray(cube[r, c])

print(f"\n  python on {SAMPLE} random scene pixels:")
print(f"  {'material':44s} {'status':>10s} {'fit>0':>7s} {'median>0':>9s} {'max':>7s}")
py = {}
for mid in group0_ids:
    result = run.evaluator.evaluate(mid, spectra)
    if result is None:
        print(f"  {mid:44s} {'None':>10s}")
        continue
    fit = np.ma.filled(np.ma.asarray(result["fit"]), 0.0).astype(float)
    py[mid] = fit
    nz = fit > 0
    med = np.median(fit[nz]) if nz.any() else 0.0
    print(f"  {mid:44s} {'ok':>10s} {int(nz.sum()):7d} {med:9.4f} {fit.max():7.4f}")

alive = sum(1 for f in py.values() if (f > 0).any())
print(f"\n  python fires on {alive}/{len(py)} evaluated group-0 materials")

print("""
Reading it:
  native fires widely + python all zero  -> group 0 is broken in Python (a
        shared preratio / udata / algorithm path, not fe3+bearing1). Chase the
        first group-0 material native maps and diff its feature trace.
  both fire                              -> group 0 is fine; the harness simply
        never compared it, and the 5821 stay unexplained -> back to bandmp.
  native silent too                      -> group 0 is genuinely dead here.
""")
#%%
import numpy as np

def check(resolved, mid, out_dir):
    basename, _ = output_info(rules[mid])
    nat = read_native_byte_gz(out_dir / f"{basename}.fit.gz").astype(int)
    fit = np.ma.filled(np.ma.asarray(resolved["fit"]), 0.0)
    both = (resolved["winner"] == mid) & (nat > 0)

    f64 = np.floor(fit.astype(np.float64) * 255 + 0.5)                   # harness today
    f32 = np.floor((fit.astype(np.float32) * np.float32(255)).astype(np.float64) + 0.5)
    frac = (fit.astype(np.float64) * 255) % 1

    for name, dn in (("float64", f64), ("float32", f32)):
        bad = both & (dn != nat)
        print(f"{mid:32s} {name}: {int(bad.sum()):4d} mismatches of {int(both.sum())}")
    bad = both & (f64 != nat)
    if bad.any():
        print(f"    fractional part of fit*255 at mismatches: "
              f"min {frac[bad].min():.6f} max {frac[bad].max():.6f}")

check(run.group_winners[4], "g4-generic.Fe2+nrw.cummingtonite", output_dirs[("group", 4)])
check(run.case_winners[1], "red.edge.shift.2", output_dirs[("case", 1)])