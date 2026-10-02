# python-tetracorder

![Cuprite, Nevada (AVIRIS 1995): 2 µm mineral map from native Tetracorder 6.00 and from python-tetracorder](https://raw.githubusercontent.com/Russjass/python-tetracorder/main/images/cuprite95_group2_comparison.jpg)
*Cuprite95 group 2 (2 µm minerals). Native Tetracorder 6.00 (centre) and python-tetracorder (right) from the same inputs,
both drawn with Tetracorder's own colouring. See [Fidelity](#fidelity).*

A port of the rules and evaluation logic of the Tetracorder algorithm to pure Python - with no guarantee or warranty!
This is merely a translation of the expertise encoded in the phenomenal
[Tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) expert system.
All citation and credit should go to:
[Roger N. Clark, Gregg A. Swayze, K. Eric Livo, Raymond F. Kokaly, Steve J. Sutley, J. Brad Dalton, Robert R. McDougal
and Carol A. Gent, 2003, *Imaging spectroscopy: Earth and planetary remote sensing with the USGS Tetracorder and expert
systems*, Journal of Geophysical Research, Vol 108.](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/documentation/clark.et.al.2003.tetracorder.JGR.2002JE001847.pdf)

*This project is an independent Python reimplementation of parts of the Tetracorder spectral identification system.
The original Tetracorder and spectroscopy software was written by Roger N. Clark and colleagues and is Copyright ©
1975-2026 Roger N. Clark, Planetary Science Institute (PSI), and contributors named in the original code. The original
software is distributed under the GNU General Public License and requires redistributions to retain its copyright
notice, conditions and disclaimer; recodings into other languages must also make the translation and source code
freely available. This project is not endorsed by Roger N. Clark, PSI, or the original contributors. The original
Tetracorder software and licensing terms are available from the PSI
[`spectroscopy-tetracorder`](https://github.com/PSI-edu/spectroscopy-tetracorder) repository.*

## Contents

- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)
- [Scene size and run time](#scene-size-and-run-time)
- [Rules, materials and reference spectra](#rules-materials-and-reference-spectra)
- [Fidelity](#fidelity)
- [Development history](#development-history)
- [Scope](#scope), [Disclaimer](#disclaimer), [License](#license)

## Installation

```bash
pip install tetracorderp
```

or, from a checkout of this repository:

```bash
pip install -e .
```

The package installs as `tetracorderp` and depends only on NumPy, SciPy, Numba, Pillow and Matplotlib. The rule set,
the reference spectra database and the post-processing lookups are bundled with it.

## Usage

`GroupEvaluator` is the entry point. Constructing it runs the full rule evaluation; the methods then write the results.

```python
from tetracorderp import GroupEvaluator

# cube: (lines, samples, bands) reflectance; wavelengths, fwhm: µm, one entry per band
run = GroupEvaluator(
    cube, wavelengths, fwhm,
    mode="default",                  # threshold preset, see tetracorderp.config.VARIABLE_PRESETS
    target_valid_bands=valid_bands,  # optional (bands,) bool, True = usable band
    masked_pixels=bad_pixels,        # optional (lines, samples) bool, True = skip this pixel
)
# Evaluation runs inside __init__: run.group_winners and run.case_winners are populated on return.

run.write_npz("scene.npz")                                    # everything, reloadable without the cube
run.write_like_tetracorder("outputs")                         # native-style .fit.gz / .depth.gz / .fd.gz
run.write_all_themes("outputs", cube_id_prefix="scene")       # colour theme maps
run.write_all_themes_display("outputs", cube_id_prefix="scene")  # theme maps with base image and key
run.write_geological_theme_images("outputs", cube_id_prefix="scene")

run = GroupEvaluator.from_npz("scene.npz")                    # later, without re-evaluating
run.export_for_qgis("qgis_outputs", cube_id_prefix="scene")
```

Other arguments:

- `temperature`, `pressure`: `(min, max)` of the scene in Kelvin and bar. Materials whose declared limits exclude the
  range are disabled, as native's `applygtpconstraints.r`.
- `disabled_groups`, `disabled_cases`, `disabled_materials`: group and case numbers, or rule ids, to skip.
- `reference_file`, `rules_file`: your own reference database or rule set, in place of the bundled ones.
- `blocking` (default `True`): evaluate in strips of about 40k pixels. Results are identical to whole-cube evaluation,
  and it is faster and uses far less memory, so there is little reason to turn it off.

### Masked pixels

`masked_pixels` follows the `numpy.ma` convention used by the EMIT and EnMAP quality flags: **True means skip**.
Masked pixels are not evaluated and come out as no winner, with zero fit and depth. Within each strip the unmasked
pixels are gathered, evaluated and scattered back, so a cloudy or partly empty scene costs only its valid pixels. The
mask is saved in the npz and restored by `from_npz`.

For example, with the EMIT L2A aggregate flag (band 7 of the mask product, 1 = flagged):

```python
masked_pixels = mask[..., 7] != 0
```

and with EnMAP L2A, combining the cloud layer with pixels that are no-data in every band:

```python
masked_pixels = (cloud != 0) | (raw == nodata).all(axis=0)
```

Flat spectra, shape `(n, bands)`, are evaluated too. The image writers need an image, so reload a flat run with
`GroupEvaluator.from_npz(path, shape=(lines, samples))`.

## Outputs

| Output | Method | Native equivalent |
|---|---|---|
| Frozen run record | `write_npz`, reloaded with `from_npz` | none |
| 8-bit fit, depth and fd per material | `write_like_tetracorder` | `<group or case dir>/<name>.{fit,depth,fd}.gz` |
| Colour theme maps | `theme_map`, `write_all_themes` | `color.results/*.png` |
| Theme maps with base image and key | `write_all_themes_display` | `color.results+labels/*+labels.png` |
| Classification and material group images | `geological_theme_images`, `write_geological_theme_images` | `geologic-classifications/`, `geologic-groups/` |
| Tetracorder for QGIS layout | `export_for_qgis` | a native run directory with geology on |

[Outputs.md](Outputs.md) traces every native product to its source, describes each python output, and lists what is
not implemented and why. The colour themes, material classes and geological lookups in `postprocessing_lookups.json`
were transcribed and in places completed by us; see [the lookups section of Outputs.md](Outputs.md#completing-the-material-components-and-other-changes-to-the-lookups).
**They are not representative of the logic contained in the authoritative Tetracorder.**

## Scene size and run time

python-tetracorder is slower than native. Indicative times on a single desktop machine:

| Scene | Pixels evaluated | Time |
|---|---:|---:|
| Cuprite95, AVIRIS (972 × 614), all enabled groups and cases | 596,808 | ~390 s |
| Drill core box scan (Specim SWIR) | 180,016 | ~120 s |
| EMIT L2A, `emit_c` mode, cloud-masked | 966,895 of 1,589,760 | 702 s |
| EnMAP L2A full cube, unmasked | 1,402,440 | 1,873 s |

Time scales with the number of pixels evaluated, so masking clouds and empty borders pays off directly.

The cube must fit in memory, as a float32 array or smaller: the evaluation converts one strip at a time and never
copies the whole cube. A scene too large for RAM, such as a full AVIRIS-NG or AVIRIS-3 flight line, should be cut into
tiles outside the package and the results stitched back together. Every pixel is evaluated independently, so tiling
does not change the results.

## Rules, materials and reference spectra

### Rules

The project is based on published/released [Tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) rule
definitions and attempts to reproduce their evaluation behaviour independently in Python. The most recent expert
system ruleset is
[cmd.lib.setup.t6.00a6](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/tetracorder.cmds/tetracorder6.00a.cmds/cmd.lib.setup.t6.00a6).
These rules have been parsed into JSON, with some manual editing, resulting in a
[dataset](tetracorderp/resources/tetracorder_rules_as_dict.json) that is read by the Python implementation.

### Reference spectra

Tetracorder rules reference spectra from USGS spectral libraries. The most recent version of the USGS spectral library
is available [here](https://www.usgs.gov/data/usgs-spectral-library-version-7-data). However, the Tetracorder rules
refer to specific reference spectra from a previous version of this library, and mapping between them is not
straightforward.

The SQLite reference spectra database in this package is derived from the libraries in the
[Tetracorder repo](https://github.com/PSI-edu/spectroscopy-tetracorder). It holds the native-resolution libraries, not
the convolved variants, with the bandpass of every sample, and the spectra are convolved to the sensor's bands at run
time using the convolver vendored from [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite).

## Fidelity

### Tetracorder 6.00 Cuprite95 validation tests

#### Native Tetracorder run

The reference run for these tests was performed in a [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite)
docker container. Tetracorder-lite vendors the entire
[authoritative Tetracorder and Specpr](https://github.com/PSI-edu/spectroscopy-tetracorder) ratfor build, with marginal
difference to the original cmd files (a renamed material, and some different comments).

The native run consumes the [Cuprite95 AVIRIS](https://popo.jpl.nasa.gov/1995_cuprite_RTGC_rfl_cube/) cube using the
Tetracorder 6.00 command/rule configuration and writes separate output products for the enabled material groups and
cases. The setup, physical conditions, disabled groups etc. are all faithfully reproduced from the
[testrun](https://github.com/PSI-edu/spectroscopy-tetracorder/tree/main/cuprite95) in the original Tetracorder, except
that this run used expert system 6.0 and the hosted run results use 5.26e1.

The following native Tetracorder internals are all recorded:

- the Cuprite95 input image cube;
- the DN offset, scale factor and deleted-data value recorded in the native run history;
- the AVIRIS95 wavelength record from the native `s06av95a` SPECpr library;
- the native `[DELETPTS]` channel mask;
- the native temperature and pressure settings;
- any groups or cases explicitly disabled by the native run;
- the native SPECpr-convolved `s06av95a` / `r06av95a` reference spectra;
- the native group and case output images used for comparison.

#### Python runs

##### Cube preparation

The Python test runs were performed on the same cube, using the same preparation as the Tetracorder internals:

- the native integer DN values are converted to float32 with exactly the native Tetracorder scaling
  `cube = (dn.astype(np.float32) + np.float32(offset)) * np.float32(scale)`;
- pixels whose raw DN equals the native deleted-data sentinel are replaced by NaN, `cube[dn == deleted_dn] = np.nan`;
- the native AVIRIS95 wavelength vector is read from
  [s06av95a](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/sl1/usgs/library06.conv/s06av95a) record 6,
  and the native `[DELETPTS]` mask is recreated.

This ensures that python-tetracorder is supplied with the same cube, wavelengths, valid-band mask, physical conditions
and rules used by the native run.

##### Library

The first test, establishing algorithmic fidelity, passed in the native convolved reference spectra. It is documented
[here](tests/fidelity_test_specpr_preconvolved.md).

The second test used the convolver vendored from [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite) on the
SQLite database derived from the Tetracorder repository. As the vendored convolver requires FWHM data per reference,
the database was updated to include that information in the samples table, and the Python code adjusted to handle the
new database format. This test established the fidelity of the vendored convolver and the internal SQLite library.
Its results are documented [here](tests/fidelity_test_full_python.md).

All further tests use the bundled database and rules, through the test script
[here](tests/Full_python_test_run_script.py). It needs the native run and the Tetracorder repository, so it cannot be
run from the package alone. The latest report, for this release, is
[tests/fidelity_test_refactor-20261002-2.md](tests/fidelity_test_refactor-20261002-2.md).

### Comparison

Native Tetracorder writes separate 8-bit products for each material, including:

- `<material>.fit.gz`
- `<material>.depth.gz`
- `<material>.fd.gz`

`GroupEvaluator.write_like_tetracorder()` emulates this output style, so the comparison evaluates both implementations
in the same output space, comparing native uint8 output directly.

Classification agreement is measured from pixels where the fit image is non-zero, while fit and depth agreement is
assessed directly from the corresponding 8-bit values.

### Results

Summary statistics are calculated for each group and case. The complete comparison output is in
[this file](tests/fidelity_test_full_python.md); only the summary is presented here.

| Group | Materials | Native only | Python only | Native pixels | Python pixels | Jaccard | Fit exact |
|---|---:|---:|---:|---:|---:|---:|---:|
| case.ep-cal-chl | 13 | 0 | 0 | 45,759 | 45,759 | 1.0000 | 1.0000 |
| case.red-edge | 2 | 0 | 0 | 528,969 | 528,969 | 1.0000 | 0.9996 |
| case.veg.type | 15 | 0 | 0 | 466,554 | 466,554 | 1.0000 | 1.0000 |
| group.1.5um-broad | 48 | 0 | 1 | 401,157 | 401,157 | 1.0000 | 1.0000 |
| group.1um | 140 | 0 | 0 | 573,079 | 573,079 | 1.0000 | 0.9999 |
| group.2um | 227 | 0 | 1 | 596,794 | 596,794 | 1.0000 | 1.0000 |
| group.2um-broad | 35 | 0 | 0 | 7,265 | 7,265 | 1.0000 | 1.0000 |
| group.ch4gas-2.3um | 29 | 0 | 0 | 6,178 | 6,178 | 1.0000 | 1.0000 |
| group.co2gas-2um | 29 | 0 | 0 | 107,050 | 107,050 | 1.0000 | 0.9997 |
| group.ree | 41 | 0 | 0 | 17,727 | 17,727 | 1.0000 | 0.9987 |
| group.ree_g21 | 40 | 0 | 0 | 17,726 | 17,726 | 1.0000 | 0.9987 |
| group.ree_samar | 29 | 0 | 0 | 6,178 | 6,178 | 1.0000 | 1.0000 |
| group.veg | 7 | 0 | 0 | 1,958,033 | 1,958,033 | 1.0000 | 0.9999 |

While there is minor degradation from the algorithmic fidelity test, incorporating the SQLite database and the
tetracorder-lite convolver has had no impact on the material classification per pixel. The differences in fit exact
(and depth exact) are attributed to floating point rounding and order of accumulation variations in the codebase,
exacerbated by the uint8 quantisation.

Native disables materials at setup and creates no output files for them (`creatoutfiles.r:76`), listing them in
`AAA.info/disabled-materials.txt`. python-tetracorder reaches the same state during evaluation: a feature whose windows
have no usable channels returns invalid, takes zero weight, and rejects the material if it is weak or must-have. Two
materials - `wollastonite_hs348.3b` and `potassium_nitrate` - therefore appear as empty Python outputs with no native
counterpart.

Every change described in [Development history](#development-history), including the packaging of this release, was
checked against this test and reproduced the table above exactly.

## Development history

### Performance pass 1 - 900+ s to ~390 s

Timings are for the full Cuprite95 fidelity run on the same machine. Each step reproduced the fidelity summary table
exactly.

| Step | Run time |
|---|---:|
| Before this pass | 900+ s |
| Caching and clean-up (@kevinhaoaus) | 781 s |
| Numba kernel for `_window_mean` | 723 s |
| `numpy.ma` removed from the evaluation path | 466 s |
| Evaluating the full cube in blocks | 382 s |
| Release clean-up: float32 strips, masked pixels, output writers | 389-397 s |

The last row is unchanged within run-to-run variation: that work was about memory, masking and the outputs, not speed.

**Caching, packaging and clean-up** - contributed by [@kevinhaoaus](https://github.com/kevinhaoaus) in
his [fork](https://github.com/kevinhaoaus/python-tetracorder). Thanks Kevin!

- **Packaging fix**: `pyproject.toml` declared `packages = ["TetracorderP"]`, which didn't match the lowercase
  `tetracorderp/` package directory. That goes unnoticed on case-insensitive filesystems such as Windows, but the
  install fails on Linux and macOS. Also added the previously undeclared `matplotlib` dependency.
- **Bug fixes**: removed two leftover debug `print()` calls in `material_evaluation.py`, and fixed a broken
  `from src.config import ...` in `tests/test_package_native_inputs_native_refs.py`.
- **Caching** (no numerical change):
  - `MaterialEvaluator.__init__` called `_prepare_rules()` twice; the duplicate call is removed.
  - `_native_not_source()` is memoized per `(target_spectra, source_rule_id, source_feature_id)`, so materials that
    veto against the same source material no longer re-evaluate it.
  - `_wtochbin()`, a pure window-to-channel lookup, is memoized by `(wavelengths, w1, w2)`.
- He also re-benchmarked the vectorised `_window_mean` that was already commented out in `tetracorder_ops.py`, and
  confirmed it is slower than the channel-by-channel loop for realistic window sizes.
- **Docs**: added the first Usage section.

**Numba kernel for `_window_mean`**

- The continuum window means are now computed by a small `@njit` kernel that sums each pixel's channels in channel
  order in REAL*4, as the Fortran DO loop does, so results are bit-identical to the previous loop.
  `error_model="numpy"` keeps an empty window giving NaN (deleted) rather than raising `ZeroDivisionError`.
- The gain was small (about 7%): the old loop was already vectorised across pixels, so most of the time was in the
  full-size temporary arrays around it rather than the loop itself.

**`numpy.ma` removed from the evaluation path**

Masked arrays allocate a data array and a mask array for every operation, and their round trips (`masked_invalid`,
`filled`, `masked_where`) copied the full pixels x channels arrays several times per feature.

- A deleted pixel is now represented by `NaN` throughout `tetracorder_ops.py` and `material_evaluation.py`, and
  `np.ma.filled(x, 0.0)` is replaced by an explicit NaN-to-zero step at the same points. Unmasked values were always
  finite (the tp1mat fit and continuum checks require it), so NaN marks exactly the pixels the mask did.
- The nvres continuum-window means reproduce `MaskedArray.mean`'s arithmetic exactly (float32 sum, then float64
  division by the count), so the red-edge case is unchanged.
- Removed six continuum-shape and band-bottom calculations from `characterise_feature` whose results were never
  read; `resolve_feature_tests` recomputes these from the current depth.
- One behavioural note: for a `relative` NOT whose relative-to feature has been deleted at a pixel, the denominator now
  takes the 1e-13 floor, matching native tp1mat (which zeroes a deleted feature). Previously it used a stale
  pre-deletion depth. This changed no pixels on Cuprite95.

**Cube evaluated in blocks**

`GroupEvaluator` evaluates the cube in strips of whole lines (`blocking=True` by default) and stitches the results back
together along the line axis.

- Every step of the pipeline is per pixel, so each strip gives exactly that strip's final answers and the stitched
  result is identical to a whole-cube evaluation. `evaluate()` takes the spectra as an argument and is run once per
  strip; `blocking=False` evaluates the whole cube in one call.
- Groups appear in every strip. A case only appears in strips where some pixel's group winner triggered it, so
  strips without it are filled with nothing detected when the results are stitched.
- The gain (466 s to 382 s, about 18%) comes from memory behaviour, not arithmetic: the per-feature pixels x channels
  temporaries and the per-material results held until group resolution are strip-sized, small enough to stay in CPU
  cache, instead of cube-sized.
- Each strip repeats the static per-feature work (window checks, reference continuum, feature weights), which offsets
  part of the gain.

The NOT-source cache was reworked to be safe with more than one set of spectra:

- It was keyed on `id(target_spectra)`, and Python reuses the ids of freed arrays, so a later strip could silently
  receive an earlier strip's NOT results.
- `MaterialEvaluator.evaluate()` now takes `same_target=False` by default and clears the cache first, so direct calls
  are always safe. `GroupEvaluator` passes `same_target=True` within a strip and calls `cache_clear()` between strips.
  The key is now just `(source_rule_id, source_feature_id)`.

The full fidelity report for this pass is in [tests/fidelity_test_blocking.md](tests/fidelity_test_blocking.md).

**Release clean-up**

- **Strip size**: strips were 64 lines; they are now about 40k pixels (`40000 // samples` lines), so wide scenes such
  as EMIT and EnMAP get strips of a similar memory size to Cuprite's.
- **float32 strips**: each strip is converted once to a contiguous float32 array, so a float64 or integer-scaled cube
  is never copied whole, and the evaluation always runs in float32 as native does.
- **Masked pixels**: `masked_pixels` added; see [Masked pixels](#masked-pixels). On EMIT it gives output identical to
  evaluating the valid pixels as a flat array and scattering them back by hand.
- **Outputs**: the theme maps, display composites with rendered keys, geological theme images and QGIS export were
  added, built on the cached 8-bit planes of `scale_like_tetracorder`; see [Outputs](#outputs).
- **Packaging**: the rules, reference database and lookups now ship inside the package as `tetracorderp/resources/`,
  found through `importlib.resources`; `GroupEvaluator` is exported from the top level.

### After this release

- Performance: evaluating strips in parallel threads is the most promising next step. Fusing more of the linear
  feature fit into Numba is estimated to give only a few per cent, as it would replace NumPy's own compiled loops.
- The geological origins images: origin vectors for the materials we added are still to be filled in.

## Scope

This repository is **not** currently intended to be:

- a drop-in replacement for the complete Tetracorder software package;
- a reproduction of its graphical tools;
- a validated operational mineral-mapping product;
- an authoritative implementation of every historical Tetracorder behaviour;
- a substitute for expert spectral interpretation.

It is an independent Python reimplementation of the parts of the algorithm required to understand, test, and reproduce
the published rule-based spectral evaluation workflow.

## Disclaimer

This software is provided for research and development purposes.

There is **no guarantee or warranty** that the implementation is:

- scientifically equivalent to the original Tetracorder implementation,
- complete,
- correct,
- or suitable for any particular application.

Do not use results from this repository as the sole basis for scientific, commercial, exploration, environmental,
safety, or other consequential decisions.

## License

python-tetracorder is licensed under the GNU General Public License, version 3 only (`GPL-3.0-only`); see
[LICENSE](LICENSE).

Third-party components keep their own licenses (details in [NOTICE](NOTICE)):

| Component | Where | License |
|---|---|---|
| Spectral convolution, from [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite) | `tetracorderp/convolve.py` | Apache-2.0 (Caltech/JPL), modified; see `licenses/LICENSE-APACHE-2.0` |
| Tetracorder command files and algorithm, from [spectroscopy-tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) | rule source, `tetracorderp/resources/` | GPL-3.0 + PSI conditions |
| USGS splib06 / sprlb06 reference spectra | `tetracorderp/resources/` | as distributed with Tetracorder |
