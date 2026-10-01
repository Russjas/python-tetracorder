"""
Top level module to evaluate sample spectra against the full Tetracorder numerical pipeline
"""
from pathlib import Path
import gzip
import re
import json
import functools
from types import SimpleNamespace

import numpy as np

from .material_evaluation import MaterialEvaluator 
from .config import NOGROUP0, VARIABLE_PRESETS, OUTPUT_DIRS
from .tetracorder_ops import prepare_rules

def to_dn(values, scale):
    """Tetracorder's 8-bit output: nint(value * scale) clipped to 0-255 (float32 product, as the Ratfor)."""
    x = np.asarray(values, dtype=np.float32) * np.asarray(scale, dtype=np.float32)
    return np.clip(np.floor(np.nan_to_num(x).astype(np.float64) + 0.5), 0, 255).astype(np.uint8)


def fd_stretch(dn):
    """Gamma stretch native applies to fd DN before display (davinci.image.to.gif -gamma), 0-255 float."""
    d = np.asarray(dn, dtype=np.float64)
    return np.clip(d * 7.5 * 0.08 ** np.sqrt(d / 400.0), 0, 255)

def _theme_rules(spec):
    """Every rule id a theme spec reads, in any of its method layouts."""
    rules = [rid for key in ("classes", "buffering_classes", "generating_classes")
             for cls in spec.get(key, []) for rid in cls["rules"]]
    rules += [channel["rule"] for channel in spec.get("channels", {}).values()]
    rules += [layer["rule"] for layer in spec.get("layers", [])]
    return rules

class GroupEvaluator:
    def __init__(self, target_spectra, target_wavelengths, target_fwhm, mode = "default", target_valid_bands = None, reference_file = None, rules_file = None,
                 disabled_groups = None, disabled_cases = None, disabled_materials = None,
                 temperature = None, pressure = None, blocking = True ):
        self.temperature = temperature    # (min, max) Kelvin of the data set, or None
        self.pressure = pressure          # (min, max) bar of the data set, or None
        self.target =  target_spectra
        self.target_shape = tuple(np.shape(target_spectra)[:-1])   # (lines, samples) image, or (n,) flat spectra
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
        if blocking:
            self.group_winners, self.case_winners = self._evaluate_by_block()
        else:
            self.group_winners, self.case_winners = self.evaluate(self.target)

    @classmethod
    def from_npz(cls, path, rules_file=None, shape=None):
        """
        Rebuild an evaluator from write_npz output, without re-running the evaluation.

        Supports the writers and theme maps (write_like_tetracorder, theme_map, write_theme_image, write_jpgs, ...).
        It holds no spectra and no reference library, so it cannot evaluate. rules_file defaults to the one recorded
        in the npz if it still exists, else the packaged rules. shape is (lines, samples) for npz files written from a
        flattened cube before write_npz recorded it.
        """
        with np.load(Path(path), allow_pickle=False) as z:
            arrays = {key: z[key] for key in z.files}
        meta = json.loads(str(arrays["meta"]))

        self = cls.__new__(cls)
        self.mode = meta["mode"]
        self.temperature = tuple(meta["temperature_K"]) if meta["temperature_K"] else None
        self.pressure = tuple(meta["pressure_bar"]) if meta["pressure_bar"] else None
        self.disabled_groups = set(meta["disabled_groups"])
        self.disabled_cases = set(meta["disabled_cases"])
        self.disabled_materials = set(meta["disabled_materials"])
        self.wavelengths = arrays["wavelengths"]
        self.valid_bands = arrays["valid_bands"]
        self.fwhm = None
        self.reference_file = meta["reference_file"]
        
        if rules_file is None and "rules_json" in arrays:
            self.rules_file = meta["rules_file"]  # provenance only; rules come from the npz
            raw_rules = json.loads(arrays["rules_json"].tobytes().decode("utf-8"))
        else:
            self.rules_file = rules_file or (meta["rules_file"] if Path(meta["rules_file"]).is_file()
                                             else self._find_file("tetracorder_rules_as_dict.json"))
            with open(self.rules_file, encoding="utf-8") as f:
                raw_rules = json.load(f)
        rules = prepare_rules(raw_rules, mode=self.mode)["rules"]
        self.evaluator = SimpleNamespace(rules=rules, disabled_materials=self.disabled_materials)

        # shape: as recorded, or as requested if it holds the same pixels
        stored = tuple(meta.get("shape") or arrays["group_1_fit"].shape)  # older npz files have no shape
        shape = tuple(shape) if shape is not None else stored
        if np.prod(shape) != np.prod(stored):
            raise ValueError(f"shape {shape} does not hold the {int(np.prod(stored))} pixels in {path}")
        self.target_shape = shape

        materials = np.array([m or None for m in arrays["materials"]], dtype=object)    # code 0 = "" = no winner
        self.group_winners, self.case_winners = {}, {}
        for key in arrays:
            match = re.fullmatch(r"(group|case)_(\d+)_material", key)
            if not match:
                continue
            kind, number = match.group(1), int(match.group(2))
            result = {"winner": materials[arrays[key]].reshape(shape)}
            for plane in ("fit", "depth", "fit_depth"):
                result[plane] = arrays[f"{kind}_{number}_{plane}"].reshape(shape)
            (self.group_winners if kind == "group" else self.case_winners)[number] = result

        self.target = np.empty(shape + (0,), dtype=np.float32)          # shape only; nothing reads the cube
        return self



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

    def evaluate(self, spectra):

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
            #print(f"Evaluating {count} of {len(self.evaluator.rules)}: {mid}")
            result = self.evaluator.evaluate(mid, spectra, same_target=True)
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
            cases.setdefault(case_num, {})[mid] = self.evaluator.evaluate(mid, spectra, same_target=True)

        #evaluate cases to find case winners
        case_winners = {}
        
        for case_num, case_results in cases.items():
            case_winners[case_num] = resolve_case(case_results, case_masks[case_num])
        return group_winners, case_winners


    def _evaluate_by_block(self, block_lines=128):
        """
        Evaluate the cube in strips of whole lines and stitch the results back together along the line axis.

        Correctness
        -----------
        Every step of the pipeline is per pixel. Evaluating a strip therefore gives exactly that strip's final answers, 
        and stacking the strips back in order gives results identical to a single whole-cube evaluate(). 

        Why this is expected to be substantially faster
        -----------------------------------------------
        The arithmetic is unchanged; the gain comes almost entirely from memory behaviour.

        1. Working arrays stop streaming from main memory.
           Each feature fit builds several pixels x channels temporaries: the valid-band-masked span, the continuum,
           the continuum-removed spectrum, and the float64 masks and products for the bandmp least-squares sums.
           Over the whole cube (i.e ~600k pixels, ~25 channels per span) each of these is tens to over a
           hundred MB, far larger than any CPU cache. Every pass over them is limited by DRAM bandwidth, and each is
           made and discarded thousands of times per run. For a strip of 64 lines (~39k pixels) the same arrays are
           a few MB, small enough to stay in L2/L3 cache between one operation and the next.

        2. Allocation gets cheaper.
           Arrays of hundreds of MB are allocated and returned to the operating system individually, and every
           fresh page is zeroed and faulted in on first touch. Strip-sized arrays are served and reused by the
           allocator without that round trip.

        3. The per-material results held until group resolution shrink.
           evaluate() keeps every group material's fit, depth, fit_depth (float64) and rejected mask until the group
           winners are resolved, then stacks them per group, which is another copy. Over the whole cube that is
           several GB of live data; per strip it is a small fraction of that, so resolution works in cache as well.
           (This is data that may later need to be retained, exposed or written TODO: consider writing intermediaries)

        4. Early exits fire more often.
           Materials skip work once no pixel is alive (e.g. the NOT vetoes stop as soon as every pixel has been
           rejected), and cases are only evaluated when some pixel triggered them. Over a whole cube there is
           almost always some pixel that keeps the work going; within a strip there often is not.

        Costs and tuning
        ----------------
        The Python-level overhead of each material and feature call (rule lookups, window checks, small-array
        bookkeeping) is paid once per strip instead of once per run, so very small strips give some of the gain
        back. block_lines trades cache fit against that repeated overhead; time a few sizes (e.g. 16, 32, 64, 128)
        on a representative cube. The only full-size arrays are the input cube (not copied when it is C-contiguous)
        and the stitched outputs, which exist twice only momentarily while each key is concatenated.

        The NOT-source cache holds pixel data for one set of spectra, and evaluate() passes same_target=True to reuse
        it across materials, so it is cleared at the start of every strip. Nothing in it could carry over anyway.
        """
        blocks = []
        for i0 in range(0, self.target.shape[0], block_lines):
            spectra = np.ascontiguousarray(self.target[i0:i0 + block_lines])    # a view when the cube is C-contiguous TODO: update usage with comment about passing contiguous in
            self.evaluator.cache_clear()   
            group_winners, case_winners = self.evaluate(spectra)
            blocks.append((spectra.shape[:-1], group_winners, case_winners))

        return (self._stitch([(shape, g) for shape, g, _ in blocks]),
                self._stitch([(shape, c) for shape, _, c in blocks]))

    @staticmethod
    def _stitch(block_results):
        """
        Stack each group's (or case's) per-strip results top to bottom. Groups appear in every strip; a case only
        appears in strips where some pixel triggered it, so missing strips are filled with nothing detected.
        """
        keys = {key for _, results in block_results for key, result in results.items() if result is not None}
        stitched = {}
        for key in keys:
            parts = {"winner": [], "fit": [], "depth": [], "fit_depth": []}
            for shape, results in block_results:
                result = results.get(key)
                if result is None:
                    result = {"winner": np.full(shape, None, dtype=object),
                              "fit": np.zeros(shape), "depth": np.zeros(shape), "fit_depth": np.zeros(shape)}
                for name in parts:
                    parts[name].append(result[name])
            stitched[key] = {name: np.concatenate(pieces, axis=0) for name, pieces in parts.items()}
        return stitched


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
    
    def _require_image(self, method):
        """Raise unless the results are an image, (lines, samples)."""
        if len(self.target_shape) != 2:
            raise ValueError(f"{method} needs image results (lines, samples), but this run has shape "
                             f"{self.target_shape}. Evaluate a (lines, samples, bands) cube, or reload with "
                             f"GroupEvaluator.from_npz(path, shape=(lines, samples)).")
    #==========writing functions =======================================================

    @functools.cached_property
    def output_specs(self):
        """rule id -> (output base name, depth/fd scale) from "8 DN 255 = x" in the rule's output block."""
        specs = {}
        for rid, rule in self.evaluator.rules.items():
            lines = [l.split("\\#")[0].strip() for l in rule["output_raw"].splitlines()]
            lines = [l for l in lines if l]
            dn, value = re.match(r"\d+\s+DN\s+(\d+)\s*=\s*(\S+)", lines[2]).groups()
            value = VARIABLE_PRESETS[self.mode].get(value, value)
            specs[rid] = (lines[1].split()[0], int(dn) / float(value))
        return specs

    def scale_like_tetracorder(self, kind, number):
        """
        One group's or case's winners as Tetracorder's 8-bit DNs, each pixel scaled by its own winner's rule:
        {"winner": object array of rule ids, "fit": uint8, "depth": uint8, "fd": uint8}, or None if no result.
        Cached; winners do not change after evaluation.
        """
        cache = self.__dict__.setdefault("_scaled", {})
        if (kind, number) not in cache:
            result = (self.group_winners if kind == "group" else self.case_winners).get(number)
            if result is None:
                cache[(kind, number)] = None
            else:
                winner = np.asarray(result["winner"], dtype=object)
                depth_scale = np.zeros(winner.shape, dtype=np.float32)
                for rid in set(winner[winner != None]):  # noqa: E711
                    depth_scale[winner == rid] = self.output_specs[rid][1]
                planes = {"fit": ("fit", 255.0), "depth": ("depth", depth_scale), "fd": ("fit_depth", depth_scale)}
                scaled = {"winner": winner}
                for plane, (key, scale) in planes.items():
                    scaled[plane] = to_dn(np.ma.filled(np.ma.asarray(result[key], dtype=np.float32), 0.0), scale)
                cache[(kind, number)] = scaled
        return cache[(kind, number)]

    @functools.cached_property
    def postprocessing_lookups(self):
        """postprocessing_lookups.json: material classes, origins and colour themes, keyed on rule id."""
        with open(self._find_file("postprocessing_lookups.json")) as f:
            return json.load(f)

    def theme_map(self, theme):
        """
        Native colour-theme map (cmds.color.support/davinci.make.*) as a (lines, samples, 3) uint8 RGB array.
        This method requires the array to have two "spatial" dimensions. This will raise if a flattened cube is passed 
        in to the GroupEvaluator constructor.

        Each rule is read as native reads its file (.fd.gif by default, .fd.gz or .depth.gz where the script does),
        0 where it did not win, from its own group or case; group-0 rules from the theme's group. The methods then
        follow their scripts: classes add a * c/255 per class, acid_buffer overlays AGM on ABM with ratio-coloured
        overlaps, rgb puts one rule per channel with a gain, binned paints full colour per DN bin.
        """
        themes = self.postprocessing_lookups["COLOUR_THEMES"]
        if theme not in themes:
            raise KeyError(f"unknown theme {theme!r}; available: {', '.join(themes)}")
        spec = themes[theme]
        self._require_image("theme_map")
        method = spec["method"]
        rules = self.evaluator.rules
        shape = self.target.shape[:2] #TODO guard here to avoid trying to write flat cubes

        # where group-0 rules are read from: the theme's one other group, else group 1 as native
        theme_rules = _theme_rules(spec)
        sources = {(rules[rid]["kind"], int(rules[rid]["number"])) for rid in theme_rules}
        others = sources - {("group", 0)}
        group0_source = others.pop() if len(others) == 1 else ("group", 1)
        sources = {group0_source if s == ("group", 0) else s for s in sources}
        if all(self.scale_like_tetracorder(kind, number) is None for kind, number in sources):
            raise ValueError(f"{theme}: no result for {sorted(sources)} (disabled, or nothing detected)")

        def image(rid, product="fd_gif"):
            """One rule's image as native reads .fd.gif / .fd.gz / .depth.gz; 0 where it did not win."""
            kind, number = rules[rid]["kind"], int(rules[rid]["number"])
            if (kind, number) == ("group", 0):
                kind, number = group0_source
            scaled = self.scale_like_tetracorder(kind, number)
            if scaled is None:
                return np.zeros(shape)                                   # native reads a missing file as 0
            if product == "fd_gif":
                plane = np.floor(fd_stretch(scaled["fd"]))               # davinci.image.to.gif -gamma, byte()
            else:
                plane = scaled[product].astype(np.float64)
            return np.where(scaled["winner"] == rid, plane, 0.0)

        def paint(classes):
            """xcolor += a * (c / 255) per class, a = sum of the class's rule images."""
            inputs = spec.get("rule_inputs", {})
            xcolor = np.zeros(shape + (3,))
            for cls in classes:
                a = sum(image(rid, inputs.get(rid, "fd_gif")) for rid in cls["rules"])
                xcolor += a[..., None] * (np.asarray(cls["rgb"], dtype=np.float64) / 255.0)
            return xcolor

        if method == "classes":
            xcolor = paint(spec["classes"])

        elif method == "acid_buffer":
            abm, agm = paint(spec["buffering_classes"]), paint(spec["generating_classes"])
            vabm, vagm = abm.sum(axis=2), agm.sum(axis=2)
            overlap = {k: np.asarray(v, dtype=np.float64) for k, v in spec["overlap"].items() if k.startswith(("r", "o"))}
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = (vagm / vabm)[..., None]
            mix = np.where(ratio < 0.33, overlap["ratio_lt_0.33"],
                        np.where(ratio > 0.66, overlap["ratio_gt_0.66"], overlap["otherwise"]))
            xcolor = abm.copy()                                          # ABM only, or nothing
            only_agm, both = (vagm > 0) & (vabm == 0), (vagm > 0) & (vabm > 0)
            xcolor[only_agm] = agm[only_agm]
            xcolor[both] = ((vagm + vabm) / 6.0)[both][:, None] * mix[both] / 255.0

        elif method == "rgb":
            xcolor = np.stack([spec["channels"][c]["gain"] * image(spec["channels"][c]["rule"]) for c in "RGB"],
                            axis=-1)

        elif method == "binned":
            xcolor = np.zeros(shape + (3,))
            for layer in spec["layers"]:
                dn = image(layer["rule"], spec.get("input", "fd_gif"))
                for b in layer["bins"]:                                  # inclusive both ends, as native
                    low, high = b["dn"]
                    xcolor[(dn > 0) & (dn >= low) & (dn <= high)] += np.asarray(b["rgb"], dtype=np.float64)

        else:
            raise NotImplementedError(f"{theme}: method {method!r}")

        return np.clip(np.floor(xcolor + 0.5), 0, 255).astype(np.uint8)  # byte(xcolor + 0.5)

    def write_theme_image(self, theme, path):
        """Write theme_map(theme) as a PNG, as native's color.results/<scene>_color-results_<theme>.png."""
        from matplotlib.image import imsave
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        imsave(path, self.theme_map(theme))

    @functools.cached_property
    def _live_groups(self):
        """Groups that write output: every group with rules, less group 0 and the disabled groups."""
        return ({int(r["number"]) for r in self.evaluator.rules.values() if r["kind"] == "group"}
                - {0} - self.disabled_groups)

    def _output_targets(self, rid):
        """(subdir, kind, number) for each directory write_like_tetracorder writes rule rid into; [] if none."""
        if rid in self.disabled_materials:
            return []
        rule = self.evaluator.rules[rid]
        number = int(rule["number"])
        if rule["kind"] == "case":
            return [] if number in self.disabled_cases else [(OUTPUT_DIRS["case"][number], "case", number)]
        if number == 0:                                                  # group 0 goes into every group that takes it
            return [(OUTPUT_DIRS["group"][g], "group", g) for g in sorted(self._live_groups - NOGROUP0)]
        return [(OUTPUT_DIRS["group"][number], "group", number)] if number in self._live_groups else []

    def write_like_tetracorder(self, output_dir):
        """
        Write every enabled material's fit, depth and fd images as Tetracorder
        does: <dir>/<name>.<plane>.gz (VICAR label + 8-bit image, gzipped) and
        <name>.<plane>.gz.hdr. Pixels a material did not win are 0. Group-0
        materials are written into every group directory that includes group 0.
        """
        output_dir = Path(output_dir)
        rules = self.evaluator.rules
        nl, ns = self.target_shape if len(self.target_shape) == 2 else (1, int(np.prod(self.target_shape)))
        lblsiz = ns if ns >= 299 else ns * (299 // ns + 1)      # creatoutfiles.r
        

        # groups that write output, and those that take group 0 (cubecorder.r:448-468)
        for rid, rule in rules.items():
            targets = self._output_targets(rid)                      # cubecorder.r:448-468
            if not targets:
                continue
            name = self.output_specs[rid][0]

            for subdir, kind, num in targets:
                (output_dir / subdir).mkdir(parents=True, exist_ok=True)
                scaled = self.scale_like_tetracorder(kind, num)

                for plane, suffix in (("fit", "FITS"), ("depth", "DEPTHS"), ("fd", "F*D")):
                    image = (np.where(scaled["winner"] == rid, scaled[plane], 0).astype(np.uint8)
                            if scaled is not None else np.zeros((nl, ns), np.uint8))

                    path = output_dir / subdir / f"{name}.{plane}"
                    title = f"{rule['output_title']:<40}{suffix}"
                    label = (
                        f"LBLSIZE={lblsiz:<16}"
                        f"FORMAT='BYTE'  TYPE='IMAGE'  BUFSIZ=20262   "
                        f"DIM=2  EOL=0  RECSIZE={ns}  ORG='BSQ'  "
                        f"NL={nl}  NS={ns}  NB=1  "
                        f"N1=0  N2=0  N3=0  N4=0  NBB=0  NLB=0  "
                        f"TASK='tetracorder'  USER='root'  "
                        f"TITLE='{title}'"
                    ).encode("ascii").ljust(lblsiz)
                    with gzip.GzipFile(f"{path}.gz", "wb", compresslevel=6, mtime=0) as f:
                        f.write(label + image.astype(np.uint8).tobytes())
                    Path(f"{path}.gz.hdr").write_text(
                        f"ENVI\ndescription = {{\n  {rule['output_title']} {suffix}\n  }}\n"
                        f"samples = {ns}\nlines   = {nl}\nbands   = 1\n"
                        f"header offset = {lblsiz}\nfile compression = 1\n"
                        "file type = ENVI Standard\ndata type = 1\ninterleave = bsq\n"
                        "sensor type = spectral data\nbyte order = 0\n"
                        "wavelength units = Micrometers\n")

    def export_for_qgis(self, output_dir, cube_id_prefix=""):
        """
        Write an output directory laid out as Grant Boxer's Tetracorder for QGIS plugin (v1.14) reads a native 6.00
        run with geology on:

            <group|case dir>/<name>.{fit,depth,fd}.gz (+ .hdr)       write_like_tetracorder
            cmds.abundances/lists.of.files.by.mineral/<list>.txt     T01 mineral layers (ABUNDANCE_LISTS)
            geologic-groups/<prefix>-<group>-group-list.txt          T02 mineral-group layers
            geologic-classifications/<prefix>-<class>-list.txt       T02 geologic-class layers

        Lists name only rules that were written. Mineral lists take each file's DN scale from its rule rather than
        native's lists, four of which are stale. Group and class lists come from RULE_LOOKUP_MATERIAL and
        MATERIAL_CLASSIFICATIONS. Georeferencing comes from the reflectance cube, chosen in the plugin.
        """
        output_dir = Path(output_dir)
        self.write_like_tetracorder(output_dir)
        lookups = self.postprocessing_lookups
        paths = {rid: [f"{subdir}/{self.output_specs[rid][0]}" for subdir, _, _ in targets]
                 for rid in self.evaluator.rules if (targets := self._output_targets(rid))}

        # T01: one list per mineral, native layout: path, DN scale, BD factor, band depth, title, library, record
        mineral_dir = output_dir / "cmds.abundances" / "lists.of.files.by.mineral"
        mineral_dir.mkdir(parents=True, exist_ok=True)
        for name, entries in lookups["ABUNDANCE_LISTS"].items():
            lines = [f"{path + '.depth.gz':<64} {255.0 / self.output_specs[e['rule']][1]:<6.4g} {e['bd_factor']:<6} "
                     f"{e['band_depth']:<7} {e['title']:<32} {e['library']} {e['record']}"
                     for e in entries for path in paths.get(e["rule"], [])]
            (mineral_dir / f"{name}.txt").write_text("".join(line + "\n" for line in lines))

        # T02: one list per mineral group and per geologic class, from the material lookups
        groups, classes = {}, {}
        lead = f"{cube_id_prefix}-" if cube_id_prefix else ""
        for rid, rule_paths in paths.items():
            for material in lookups["RULE_LOOKUP_MATERIAL"].get(rid) or []:
                info = lookups["MATERIAL_CLASSIFICATIONS"].get(material) or {}
                if info.get("material_group"):
                    groups.setdefault(info["material_group"], set()).update(rule_paths)
                if info.get("Classification"):
                    classes.setdefault(info["Classification"], set()).update(rule_paths)
        for folder, table, suffix in (("geologic-groups", groups, "-group"), ("geologic-classifications", classes, "")):
            (output_dir / folder).mkdir(exist_ok=True)
            for label, rule_paths in table.items():
                lines = sorted(rule_paths) if suffix else [f"{label:>20}   {p}" for p in sorted(rule_paths)]
                (output_dir / folder / f"{lead}{label}{suffix}-list.txt").write_text("".join(l + "\n" for l in lines))



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
            rules_json           uint8, the rules file the run used (UTF-8 bytes)
            meta                 JSON string: mode, files, conditions, disabled lists

        Read back with np.load(path); no pickling is needed or allowed.
        """
        materials = [""] + list(self.evaluator.rules)
        code = {mid: i for i, mid in enumerate(materials)}
        arrays = {"materials": np.array(materials, dtype=str),
                "wavelengths": np.asarray(self.wavelengths, dtype=np.float32),
                "valid_bands": np.asarray(self.valid_bands, dtype=bool),
                 "rules_json": np.frombuffer(Path(self.rules_file).read_bytes(), dtype=np.uint8)}

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
            "shape": list(self.target_shape),
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