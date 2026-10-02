# Outputs

## Native tetracorder products - from spectroscopy-tetracorder
Traced from `tetracorder.cmds/tetracorder6.00a.cmds/` (`cmd.runtet`, `cmds.all.support`, `cmds.color.support`,
`cmds.geology`, `cmds.abundances`) and `tetracorder6.00/cubecorder.r`.

Every product is the same pattern: **a list of rules → each rule's winner image → combine → render.** They differ only
in which list, which image (fit, depth, fd, or the gamma-stretched fd.gif), how images are combined, and the output.

| My name | Product | Made by | What it is |
|---|---|---|---|
| 8-bit data | `<group or case dir>/<name>.{fit,depth,fd}.gz` (+ `.hdr`) | the run (`cubecorder.r`) | One 8-bit VICAR image per rule per plane, non-zero only where that rule won its group or case. The material is identified by the file name. Group-0 rules are copied into every group directory that takes group 0; case rules go in case directories. DN = nint(v·255/x) for depth and fd, from the rule's `8 DN 255 = x`; fit DN = nint(fit·255). |
| fd-gamma | `<name>.fd.gif` | `cmds.all.support/gen.fd.gif.images` | fd DN with gamma stretch `a·7.5·0.08^sqrt(a/400)`, byte. Mostly an intermediate for themed images and geological themed images. |
| fd-overlays | `results.dual.<dir>/<name>.fd.ovrly.dual.gstr.jpg` | `gen.fd.jpg.overlay+base-dual.images` (`davinci.red.overlay.gray-image … -25 … -gamma -dual`) | Each rule's fd as a red overlay on the grey base band (`base-image/base-image.jpg`), autoscale limit −25, gamma, with the plain base beside it. One per rule per directory. A `-notzero` follow-up keeps only rules that won somewhere. Seems to be bugged in native |
| Themed images | `color.results/*.png`, `color.results+labels/*+labels.png` | `cmds.color.support/davinci.make.*`, `make.color.results.all`, `davinci.cat.colormap+key-below` | Theme maps (27 scripts; 21 produced on Cuprite 95). Per class: sum the rules' `.fd.gif` (a few read raw `.fd.gz` or `.depth.gz`), add `a·c/255`, write `byte(xcolor+0.5)`. Exceptions: acid-minerals (ABM/AGM overlap coloured by intensity ratio), veg-water-rgb (one rule per channel, gains 1, 1.5, 2), red-edge-shift (full colour per depth-DN bin). Base image (`base-image/color-visRGB.jpg`, from `COLOR.channels`, else `base-image.jpg`) pasted to the right; hand-drawn key PNG below for `+labels`. |
| geological origins data | `geologic-origins/geologic_origins_cube.v` + 11 `…chan-NN-<origin>.gif` | the run, when `geology` is set (`cubecorder.r:745–885`); `cmd.make.gifs.from.geocube` | Per pixel, for each enabled group's winner m and each of its material components c: `w = min(1, 1.6·|depth|)`, `x = origin[c,k]/9 · fit · confidence/9 · w`; if x > 1e-4, numerator += x and denominator += \|depth\| + 0.2. Cube = numerator/denominator, written int16 ×1000, BIL, 11 bands. Each channel → gif as byte(DN × 0.254), no gamma. Cases not included. |
| Geological theme images | `geologic-classifications/<cubeid>-<class>-list.txt` + `-class.gif`; `geologic-groups/<cubeid>-<group>-group-list.txt` + gif | `cmds.geology/cmd.make-all.*`, called by `cmd.runtet` (lines 888–900) | One list and one greyscale image per material class (or mineral group): byte(sum of `.fd.gif` of every rule with a component in that class). Lists are built from `AAA.info/material-classifications.txt`, written by the run. Class match is a substring grep. |
| Masses | `results.masses/<mineral>-{bdepth,relbdepth,mass…}` | `cmds.abundances/cmd.compute.abundances` → `davinci.image.material.mass`, called by `cmd.runtet` | 17 EMIT dust minerals. Per mineral, a hand-curated list (`lists.of.files.by.mineral/*.txt`: depth file, DN scale, BD factor, reference band depth, title, library, record); band depth = Σ DN/255·scale·factor, capped; relative depth ÷ 0.355 (kaolinite, for every mineral); mass by Clark & Roush random walk with per-mineral Δk and grain size, density 3.1 for every mineral; slab-model mass. Typo: the goethite-fine command runs `delk 5125120 0.0001` while echoing `delk 5120 gd 0.0001`. Four list entries carry a stale DN scale (maghemite, illite.roscoelite: list 0.5, rule 1.0). |
| Abundances | `results.abundances/model1–4/` | `cmd.compute-model-abundances`, by hand | EMIT-specific abundance models (linear band depth, Beer's-law slab, random-walk mass, model 4). Needs the cube and optical constants. |

## Python outputs

| My name | native equivalent | python call | Made by | What it is |
|---|---|---|---|---|
| .npz | No | `GroupEvaluator.write_npz(filepath.npz)` | numpy serialisation | A frozen record of the results of a python-tetracorder evaluation run: per group and case, the winner (as material codes) plus fit, depth and fit_depth as float32, and the wavelengths and valid bands. Also embeds the rules file the run used, the image shape and the base image, the pixel mask if one was given, and records the mode, temperature, pressure and disabled sets in `meta`. Reloaded with `GroupEvaluator.from_npz(path)`, which supports every writer below but cannot re-evaluate. |
| Base image | `base-image/color-visRGB.jpg`, `base-image/base-image.jpg` | `GroupEvaluator.base` (attribute); `GroupEvaluator.base_image(wavelengths)` | `fd_stretch` | The backdrop for the Themed images display. Built once when the run is evaluated, reading only the three chosen bands (`None` for flat runs), stored in the .npz and restored by `from_npz`. An equivalent of native's, not a reproduction: the valid bands nearest true colour (0.64/0.55/0.47 µm), else SWIR (2.20/1.65/1.25 µm), else the first, middle and last valid channels; a shared 1 % low clip, a linear stretch, then fd-gamma's tone curve. Set `run.base = run.base_image(...)` for other bands. |
| 8-bit data (in memory) | 8-bit data | `GroupEvaluator.scale_like_tetracorder(kind, number)` | rule `8 DN 255 = x` via `output_specs`; `fd_stretch` | The native DN planes for one group or case, as arrays rather than files: `winner`, `fit`, `depth`, `fd` (uint8), each scaled by its pixel's winning rule, plus `fd_gamma`, the in-memory equivalent of native's `.fd.gif`. Cached. Every image writer below is built on it. |
| 8-bit data | 8-bit data | `GroupEvaluator.write_like_tetracorder(output_dir)` | `scale_like_tetracorder`, VICAR label + gzip, ENVI `.hdr` | Native's `<group or case dir>/<name>.{fit,depth,fd}.gz` and `.gz.hdr`, including group 0 rules copied into every directory that takes group 0. Matches native: Jaccard 1.0000 in all 13 Cuprite 95 directories, fit exact ≥ 0.9987. A flat (non-image) run is written as a one-line image. |
| fd-gamma | fd-gamma | none (internal) | `scale_like_tetracorder`, `fd_stretch` | Not written as files, by decision. Held as the `fd_gamma` plane of `scale_like_tetracorder` and used wherever native reads a `.fd.gif`. |
| fd-overlays | fd-overlays | not implemented | | Not implemented, by decision. They add nothing beyond the 8-bit data and the Themed images. Native's are empty on our reference run because of the base-image bug in `gen.fd.jpg.overlay+base-dual.images`. |
| Themed images | `color.results/*.png` (colour map) | `GroupEvaluator.theme_map(theme)` → array; `write_all_themes(output_dir, cube_id_prefix="")` | `COLOUR_THEMES`; `scale_like_tetracorder` | All 27 native themes (`classes`, `acid_buffer`, `rgb`, `binned`) as colour maps alone. Each rule is read as native reads its file (`.fd.gif`, or `.fd.gz`/`.depth.gz` per `rule_inputs`). Group 0 rules come from the theme's group, or group 1 if every rule is group 0. `write_all_themes` writes `color.results/[<prefix>_]color-results_<theme>.png` and skips themes with no result. Image runs only. |
| Themed images (display) | `color.results+labels/*+labels.png` | `GroupEvaluator.write_all_themes_display(output_dir, cube_id_prefix="", base=None)` | `theme_map`, `self.base`, `postprocessing.display_legend` (Pillow, DejaVu Sans Bold) | The colour map with the base image to its right and a key below, as `color.results+labels/[<prefix>_]color-results_<theme>+labels.png`. The key is rendered from `COLOUR_THEMES` in native's style: white bold on black, wrapped title, section columns, dark-to-full gradient swatches (flat for bins), "python-tetracorder" footer. Scenes narrower than 600 px get the key as a single column to the right instead. Native's dual-only image is not produced. Image runs only. |
| Geological origins data | Geological origins data | images planned | | The 11-channel cube is not written, by decision. Only the 11 per-origin images derived from it are planned, computed in memory with native's weighting from `GEOLOGICAL_ORIGIN_CHANNELS` and `RULE_LOOKUP_MATERIAL`. Curated mapping, so not a faithful reproduction. |
| Geological theme images | Geological theme images | `GroupEvaluator.geological_theme_images()` → arrays; `GroupEvaluator.write_geological_theme_images(output_dir, cube_id_prefix="")` | `RULE_LOOKUP_MATERIAL` → `MATERIAL_CLASSIFICATIONS`; `scale_like_tetracorder`, `fd_stretch` | One greyscale image per `Classification` and per `material_group` value, in one pass: `byte(Σ fd-gamma)` over every rule_id with a material carrying that value, in every directory it is written to. Written with their lists, as native names them: `geologic-classifications/[<prefix>-]<Classification>-class.gif` and `-list.txt`, and `geologic-groups/[<prefix>-]<material_group>-group.gif` and `-group-list.txt`. Values with nothing detected get a list but no image, as native. Curated mapping. Image runs only. |
| Masses | Masses | out of scope | | EMIT-specific post-process. The per-mineral lists are kept as `ABUNDANCE_LISTS`, used by `export_for_qgis`. |
| Abundances | Abundances | out of scope | | Needs the reflectance cube and optical constants; EMIT-specific. |
| QGIS export | No (mimics a native 6.00 run directory with geology on) | `GroupEvaluator.export_for_qgis(output_dir, cube_id_prefix="")` | `write_like_tetracorder`; `ABUNDANCE_LISTS`; `_write_geological_theme_lists` | A shim for Grant Boxer's Tetracorder for QGIS plugin (v1.14). Writes a more authentic Tetracorder output directory, with the expected book-keeping files: the 8-bit data, the 21 mineral lists (DN scale from the rules) and the geological theme lists. **Plugin not yet published - citation and details to follow** |



## Outputs not implemented

python-tetracorder does not produce fd-gamma or fd-overlays files, the geological origins data, the masses or the
abundance models.

- **fd-gamma and fd-overlays, by decision.** fd-gamma is applied in memory wherever native reads a `.fd.gif`, so
  written files would add nothing. The fd-overlays add nothing beyond the 8-bit data and the Themed images. Native's
  own fd-overlays are empty on our reference run anyway, because of a bug in `gen.fd.jpg.overlay+base-dual.images`
  (the base image is written with a doubled `.base.jpg` and copied from the wrong directory).
- **The geological origins data (the 11-channel cube), by decision.** Only the 11 per-origin images derived from it
  are wanted, and those are still to be written. Their inputs (`GEOLOGICAL_ORIGIN_CHANNELS`, `RULE_LOOKUP_MATERIAL`)
  are in place, but they depend on our curated material mapping, so they will not reproduce native pixel for pixel.
- **The masses, out of scope.** They are an EMIT-specific post-process with hard-coded constants: relative depth to
  kaolinite for every mineral, density 3.1, and one mistyped Δk. Their per-mineral lists are kept as
  `ABUNDANCE_LISTS` and used by `export_for_qgis`, but the physics is not ported.
- **The model abundances, out of scope.** They need the reflectance cube and optical constants.

Two outputs are produced differently from native:

- **The base image** behind the Themed images display is an equivalent of native's, not a reproduction. It is built
  once when the run is evaluated and stored in the .npz.
- **The keys** are rendered from `COLOUR_THEMES` in native's style, rather than copied from native's hand-drawn key
  PNGs.  

## Completing the material components, and other changes to the lookups

Native assigns each rule's components from the `list materials` block of `cmd.lib.setup.t6.00a6`, then matches each
component by its first 30 characters against `material-classification.txt` and
`geological-origin-table-by-material.txt`. An unmatched component falls to the tables' `zero` row and drops out of
every class and origin. That chain is incomplete:

- **Blank rules.** 248 of the 707 rules have an empty `list materials` block, so they get no class and no origins.
  All but one (`wollastonite_hs348.3b`) come after a `ZZZZZZ edit to here` marker in the cmd file. They cover all of
  groups 11, 13, 14, 15, 19, 24, 26, 27, 37 and 38 (bar a few rules) and every case. The classification table
  already holds rows only these rules would use (antigorite, lizardite, endellite, serpentine, tremolite, uralite,
  smaragdite and others), so this reads as unfinished work rather than intent.
- **Components that match no row.** 72 component lines in 49 + 23 rules miss the tables:
  - prefix and spelling drift: `feldspar_albite`, `pyroxene_augite`, `yroxene.bronzite`, `hemtatite`, `hypersthene`
    vs `pyroxene_hypersthene`;
  - wrong entries: `phyllosilicate` written as palygorskite's component, and `ZZZZZZ edit to here` as benzene's;
  - materials with no row at all: vermiculite, magnetite, desert varnish, copper sulfate, copper precipitate,
    generic soil, and `neodymium`, used by 20 Nd-bearing minerals.
- **The two tables disagree.** `strontinite` in the origins table vs `strontianite` everywhere else; alunite has a
  class but no origins; aliphatic, manganese and propylzone have origins but no class.
- **Parse errors.** A few component lines are written as `name abundance`, which native reads as `name class`, so
  the abundance becomes the class.

What we did:

- **We left the rules untouched,** and built `RULE_LOOKUP_MATERIAL`, keyed by rule id, mapping each rule to one or
  more materials. Every one of the 707 rules now resolves to entries in `MATERIAL_CLASSIFICATIONS`. There is no
  `zero` fallback.
- **Typos and prefix drift were mapped** to the obvious existing row (e.g. `feldspar_albite` → `albite`,
  `hemtatite` → `hematite`).
- **Two new materials were added:** `vermiculite`, and `neodymium` as its own material. The Nd-bearing minerals are
  not folded into `neodymium_oxide`, which would class monazite and apatite as oxides. Other new keys were added
  where they were uncontroversial.
- **Blank rules were filled by judgement,** using the rule ids, output titles and reference spectra, and reviewed in
  CSV tables with notes. The existing serpentine and amphibole rows were used for their group 13 rules. These are
  judgement calls, not reproductions of native.
- **Known gap:** origin vectors for the newly added materials are still to be filled in
  (`GEOLOGICAL_ORIGIN_CHANNELS` has 213 materials against 246 classifications), so those materials contribute
  nothing to an origins product.

**NB** This extension of the rules components and materials list was performed exclusively by us, with all mistakes in judgement our own.  
As such these look up tables are not representative of the logic contained in the authoritative Tetracorder.

The other tables in `postprocessing_lookups.json`:

- **`COLOUR_THEMES`:** the 27 native theme scripts (`cmds.color.support/davinci.make.*`) transcribed to data:
  - method (`classes`, `acid_buffer`, `rgb`, `binned`), class colours and rules;
  - labels transcribed from the hand-drawn key images;
  - script file names mapped to rule ids;
  - `rule_inputs` where a script reads raw `.fd.gz` or `.depth.gz` instead of the stretched `.fd.gif`
    (1micron-minerals-a's Fe²⁺ class; nacrite in the three 2-micron themes).
- **`ABUNDANCE_LISTS`:** native's 21 `lists.of.files.by.mineral` files, keyed by rule id:
  - paths and DN scales are dropped and regenerated from the rules on export, which fixes four stale scales
    (maghemite, illite.roscoelite);
  - inline `# notes` are stripped, and `AAA.Readme.EMIT.txt`, which sits among them, is excluded.
- **`GEOLOGICAL_ORIGIN_NAMES` / `GEOLOGICAL_ORIGIN_CHANNELS`:** native's 11 origin channels and the per-material
  weight vectors, carried over from `geological-origin-table-by-material.txt`.


