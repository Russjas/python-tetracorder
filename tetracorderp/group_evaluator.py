"""
Top level module to evaluate sample spectra against the full Tetracorder numerical pipeline
"""
from pathlib import Path
import gzip
import re
import json
import functools

import numpy as np

from .material_evaluation import MaterialEvaluator 
from .config import NOGROUP0, VARIABLE_PRESETS, OUTPUT_DIRS

def to_dn(values, scale):
    """Tetracorder's 8-bit output: nint(value * scale) clipped to 0-255 (float32 product, as the Ratfor)."""
    x = np.asarray(values, dtype=np.float32) * np.asarray(scale, dtype=np.float32)
    return np.clip(np.floor(np.nan_to_num(x).astype(np.float64) + 0.5), 0, 255).astype(np.uint8)


def fd_stretch(dn):
    """Gamma stretch native applies to fd DN before display (davinci.image.to.gif -gamma), 0-255 float."""
    d = np.asarray(dn, dtype=np.float64)
    return np.clip(d * 7.5 * 0.08 ** np.sqrt(d / 400.0), 0, 255)



class GroupEvaluator:
    def __init__(self, target_spectra, target_wavelengths, target_fwhm, mode = "default", target_valid_bands = None, reference_file = None, rules_file = None,
                 disabled_groups = None, disabled_cases = None, disabled_materials = None,
                 temperature = None, pressure = None, blocking = True ):
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
        if blocking:
            self.group_winners, self.case_winners = self._evaluate_by_block()
        else:
            self.group_winners, self.case_winners = self.evaluate(self.target)

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

        Each pixel takes the summed RGB of every class listing its winning rule (mixtures blend, as emit8),
        scaled by the winner's gamma-stretched fd DN; rules in the theme's "depth_input" use their raw depth
        DN instead, as native reads those from .depth.gz. Pixels whose winner is not in the theme are black.

        Available themes: '1micron-minerals-a', 'hematite+goethite.grain.size-a', 'water-a', 'veg,water,snow', 
        'snow-grain-size-water.a', '2micron-minerals', '2micron-minerals-b4', '2micron-minerals-detail2', '2micron-mins-emit8', 
        '2micron-minerals-muscovite-comp', 'prehnite-chlorite-mix+perchlorate', 'organics-veg-2um-a', 'pyroxene.2um.band.position', 
        '1.5um.broadfeats', '1.9um.water.wave.position.a', '1.9um.water.band.position', '1.9um.water.sulfates.band.position', 
        '1.9um.water.zeolites.band.position', '2.8um.oh.band.position', '3um.waterfeats', '3.5um.feat.position', 'ree.b-g21', 
        'vegetation-cover-a', 'veg-spectral-type', 'acid-minerals-buffering-minerals.a', 'veg-water-rgb', 'red-edge-shift-a'
        """
        if theme not in self.postprocessing_lookups["COLOUR_THEMES"].keys():
            raise ValueError(f"{theme} is not an accepted theme")
        spec = self.postprocessing_lookups["COLOUR_THEMES"][theme]
        if spec["method"] != "classes":
            raise NotImplementedError(f"{theme}: method {spec['method']!r} not implemented yet")

        rule_rgb = {}
        for cls in spec["classes"]:
            for rid in cls["rules"]:
                rule_rgb[rid] = rule_rgb.get(rid, 0) + np.asarray(cls["rgb"], dtype=np.float64)

        # the theme reads one group or case; group-0 rules arrive as winners inside it
        rules = self.evaluator.rules
        
        sources = {(rules[rid]["kind"], int(rules[rid]["number"])) for rid in rule_rgb}
        if len(sources) > 1:
            sources -= {("group", 0)}             # group 0 rules win inside whichever group the theme is bounded to
        if sources == {("group", 0)}:
            sources = {("group", 1)}              # native reads group-0-only themes (water-a, snow...) from group.1um
        if len(sources) != 1:
            raise ValueError(f"{theme}: rules span {sorted(sources)}, expected one group or case")
        
        kind, number = sources.pop()
        scaled = self.scale_like_tetracorder(kind, number)
        if scaled is None:
            raise ValueError(f"{theme}: {kind} {number} has no result (disabled, or nothing detected)")

        winner = scaled["winner"]
        strength = np.floor(fd_stretch(scaled["fd"]))                 # davinci byte() of the .fd.gif
        for rid in spec.get("depth_input", []):
            won = winner == rid
            strength[won] = scaled["depth"][won]                      # native reads these from .depth.gz

        rgb = np.zeros(winner.shape + (3,), dtype=np.float64)
        for rid, colour in rule_rgb.items():
            won = winner == rid
            if won.any():
                rgb[won] = strength[won][:, None] * colour / 255.0
        return np.clip(np.floor(rgb + 0.5), 0, 255).astype(np.uint8)  # byte(xcolor + 0.5)

    def write_theme_image(self, theme, path):
        """Write theme_map(theme) as a PNG, as native's color.results/<scene>_color-results_<theme>.png."""
        from matplotlib.image import imsave
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        imsave(path, self.theme_map(theme))


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
                targets = [(OUTPUT_DIRS["case"][number], "case", number)]
            elif number == 0:
                targets = [(OUTPUT_DIRS["group"][g], "group", g) for g in takes_group0]
            elif number in live_groups:
                targets = [(OUTPUT_DIRS["group"][number], "group", number)]
            else:
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


    def write_jpgs(self, output_dir, dpi=150):
        from matplotlib.figure import Figure
        from matplotlib.patches import Patch
        from matplotlib import colormaps
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        rules = self.evaluator.rules
        nl, ns = self.target.shape[:2]
        cube = np.asarray(self.target, dtype=np.float32)[..., np.asarray(self.valid_bands, dtype=bool)]
        grey = np.nanmean(np.where(cube > -1e30, cube, np.nan), axis=2)
        lo, hi = np.nanpercentile(grey, (2, 98))
        grey = np.nan_to_num(np.clip((grey - lo) / (hi - lo), 0, 1))
        palette = np.vstack([colormaps[c].colors for c in ("tab20", "tab20b", "tab20c")])
        palette = palette[np.ptp(palette, axis=1) > 0.1]          # greys would vanish into the background
        for kind, results in (("group", self.group_winners), ("case", self.case_winners)):
            for number, result in sorted(results.items()):
                if result is None:
                    continue
                winner = np.asarray(result["winner"], dtype=object)
                ids, counts = np.unique(winner[winner != None].astype(str), return_counts=True)  # noqa: E711
                if ids.size == 0:
                    continue
                order = np.argsort(counts)[::-1]
                rgb = np.repeat(grey[..., None], 3, axis=2)
                handles = []
                slots = [m for m, r in rules.items() if r["kind"] == kind
                        and int(r["number"]) in ({number, 0} if kind == "group" else {number})]
                fd_dn = self.scale_like_tetracorder(kind, number)["fd"]
                for mid, n in zip(ids[order], counts[order]):
                    colour = palette[slots.index(mid) % len(palette)]
                    won = winner == mid
                    name = self.output_specs[mid][0]
                    w = (fd_stretch(fd_dn[won]) / 255.0)[:, None]
                    rgb[won] = rgb[won] * (1 - w) + colour * w      # gen.fd.gif.images gamma
                    handles.append(Patch(color=colour, label=f"{name}  ({n} px)"))
                fig = Figure(figsize=(8 * ns / nl + 3, 8))
                ax = fig.subplots()
                ax.imshow(np.clip(rgb, 0, 1), interpolation="nearest")
                ax.set_axis_off()
                ax.set_title(f"{kind} {number}: {OUTPUT_DIRS[kind][number]}")
                ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                            fontsize=7, frameon=False)
                fig.savefig(output_dir / f"{kind}{number:02d}_{OUTPUT_DIRS[kind][number].split('.', 1)[1]}.jpg",
                            dpi=dpi, bbox_inches="tight", pil_kwargs={"quality": 95})

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
                   fd_dn = self.scale_like_tetracorder(kind, number)["fd"]
                   rgb = np.zeros((nl, ns, 3))
                   for mid in set(winner[winner != None]):  # noqa: E711
                       won = winner == mid
                       w = (fd_stretch(fd_dn[won]) / 255.0)[:, None]
                       rgb[won] = palette[slots.index(mid) % len(palette)] * w      # gen.fd.gif.images gamma
                   imsave(output_dir / f"{kind}{number:02d}_{OUTPUT_DIRS[kind][number].split('.', 1)[1]}.jpg",
                          np.clip(rgb, 0, 1), pil_kwargs={"quality": 95})

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