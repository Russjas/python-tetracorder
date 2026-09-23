# python-tetracorder

![Cuprite, Nevada (AVIRIS 1995): 2 µm mineral map from native Tetracorder 6.00 and from python-tetracorder](resources/cuprite95_group2_comparison.jpg)
*Cuprite95 group 2 (2 µm minerals). Native Tetracorder 6.00 (centre) and python-tetracorder (right) from the same inputs, both drawn with Tetracorder's own colouring. See [Fidelity](#fidelity).*

A port of the rules and evaluation logic of the tetracorder algorithm to pure Python - with no guarantee or warranty!  
This is merely a translation of the expertise encoded in the phenomenal [Tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) expert system.  
All citation and credit should go to:  
[Roger N. Clark, Gregg A. Swayze, K. Eric Livo, Raymond F. Kokaly, Steve J. Sutley, J. Brad Dalton, Robert R. McDougal and Carol A. Gent, 2003, *Imaging spectroscopy: Earth and planetary remote sensing with the USGS Tetracorder and expert systems*, Journal of Geophysical Research, Vol 108.](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/documentation/clark.et.al.2003.tetracorder.JGR.2002JE001847.pdf)

*This project is an independent Python reimplementation of parts of the Tetracorder spectral identification system. The original Tetracorder and spectroscopy software was written by Roger N. Clark and colleagues and is Copyright © 1975-2022 Roger N. Clark, Planetary Science Institute (PSI), and contributors named in the original code. The original software is distributed under the GNU General Public License and requires redistributions to retain its copyright notice, conditions and disclaimer; recodings into other languages must also make the translation and source code freely available. This project is not endorsed by Roger N. Clark, PSI, or the original contributors. The original Tetracorder software and licensing terms are available from the PSI [`spectroscopy-tetracorder`](https://github.com/PSI-edu/spectroscopy-tetracorder) repository.*


## Goals

The project is intended to:

- port Tetracorder mineral identification rules into a structured, machine-readable form;
- reproduce the relevant spectral evaluation logic in pure Python;
- make the rule evaluation process easier to inspect, test, modify, and integrate into other Python workflows;
- provide a foundation for testing Tetracorder-style mineral identification against modern hyperspectral datasets;
- The implementation is therefore intended to rely primarily on standard Python scientific tools such as:
  - NumPy
  - SciPy
  - SQLite

## Rules/Materials

The project is based on published/released [Tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) rule definitions and attempts to reproduce their evaluation behaviour independently in Python.
The most recent expert system ruleset is [cmd.lib.setup.t6.00a6](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/tetracorder.cmds/tetracorder6.00a.cmds/cmd.lib.setup.t6.00a6). These rules have been parsed into json format, with some manual editing resulting in a [dataset](resources/tetracorder_rules_as_dict.json) that can be read by the python implementation  

## Reference spectra

Tetracorder rules reference spectra from USGS spectral libraries. The most recent version of the USGS spectral library is available [here](https://www.usgs.gov/data/usgs-spectral-library-version-7-data). However, the Tetracorder rules refer to specific reference spectra from a previous version of this library, and mapping between them is not straightforward.
The SQLite reference spectra database in this repository is derived from the libraries in the [Tetracorder repo](https://github.com/PSI-edu/spectroscopy-tetracorder). The database here is derived from the native libraries rather than the convolved variants.

## Fidelity

### Tetracorder 6.00 Cuprite95 validation tests

#### Native Tetracorder run

The reference run for these tests was perfomed in a [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite) docker container. Tetracorder-lite vendors the entire [authoritative Tetracorder and Specpr](https://github.com/PSI-edu/spectroscopy-tetracorder) ratfor build, with marginal difference to the original cmd files (a renamed material, and some different comments).

The native run consumes the [Cuprite95 AVIRIS](https://popo.jpl.nasa.gov/1995_cuprite_RTGC_rfl_cube/) cube using the Tetracorder 6.00 command/rule configuration and writes separate output products for the enabled material groups and cases. The setup, physical conditions, disabled groups etc are all faithfully reproduced from the [testrun](https://github.com/PSI-edu/spectroscopy-tetracorder/tree/main/cuprite95) in the original Tetracorder, except this run used expert system 6.0 and the hosted run results use 5.26e1.

Native Tetracoder internals:

- the Cuprite95 input image cube;
- the DN offset, scale factor and deleted-data value recorded in the native run history;
- the AVIRIS95 wavelength record from the native `s06av95a` SPECpr library;
- the native `[DELETPTS]` channel mask;
- the native temperature and pressure settings;
- any groups or cases explicitly disabled by the native run;
- the native SPECpr-convolved `s06av95a` / `r06av95a` reference spectra;
- the native group and case output images used for comparison.

are all recorded.


#### Python runs

##### Cube preparation

The python side test runs were perfomed on the same cube, using the same preparation as the Tetracorder internals.
 - converts the native integer DN values to float32 and applies exactly the native Tetracorder scaling  
 ```cube = (dn.astype(np.float32) + np.float32(offset)) * np.float32(scale)```  
 - Pixels whose raw DN equals the native deleted-data sentinel are replaced by NaN  
```cube[dn == deleted_dn] = np.nan```
 - The resulting cube is cached as a .npy, then loaded and made contiguous  
 ```cube = np.ascontiguousarray(full_cube[LINES])```
  - reads the native AVIRIS95 wavelength vector from [s06av95a](https://github.com/PSI-edu/spectroscopy-tetracorder/blob/main/sl1/usgs/library06.conv/s06av95a) record 6 and recreates the native [DELETPTS] mask.

This ensures that `python-tetracorder` is supplied with the same cube, wavelengths, valid-band mask, physical conditions and rules used by the native run.

##### Library

The first test to establish algorithmic fidelity was ran by passing in the native convolved reference spectra. It is documented [here](tests/fidelity_test_specpr_preconvolved.md)  

The second test used a vendored convolve function from [tetracorder-lite](https://github.com/emit-sds/tetracorder-lite) on the SQLite database derived from the Tetracorder repository.  

As the vendored convolver requires fwhm data per reference, the database was updated to include that information in the samples table. The python code has also been adjusted to handle the new database format.

All python internals were used in this test; the SQLite database, the python convolver, `GroupEvaluator`, `MaterialEvaluator` and the `write_like_tetracorder()` method.

The full test script is documented [here](tests/Full_python_test_run_script.py)

All enabled groups and cases are evaluated. Group-0 materials are included in groups where the native configuration includes group 0.

### Comparison

Native Tetracorder writes separate 8-bit products for each material, including:

 - `<material>.fit.gz`
 - `<material>.depth.gz`
 - `<material>.fd.gz`

The `write_like_tetracorder()` method on `GroupEvaluator` emulates this output style.
The comparison therefore evaluates both implementations in the same output space; comparing native uint8 output directly.

Classification agreement is measured from pixels where the fit image is non-zero, while fit and depth agreement is assessed directly from the corresponding 8-bit values.

### Results  

The python implementation is significantly slower - although there are still optimisations that can be performed once behaviour is equivalent.

Summary statistics are also calculated for each group and case. 

The complete test comparison output is in [this file](tests/fidelity_test_full_python.md). Only a summary is presented here.


#### Summary


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

While there is minor degradation from the algorithmic fidelity test, incorporating the SQLite database and the tetrapy convolver, has had no impact on the material classification per pixel. The difference in fit exact (and depth exact) are attributed to floating point rounding, order of accumulation variations in the codebase, and exacerbated by the uint8 quantisation.  

Native disables materials at setup and creates no output files for them (creatoutfiles.r:76), listing them in AAA.info/disabled-materials.txt. python-tetracorder reaches the same state during evaluation: a feature whose windows have no usable channels returns invalid, takes zero weight, and rejects the material if it is weak or must-have. Two materials - wollastonite_hs348.3b and potassium_nitrate - therefore appear as empty python outputs with no native counterpart. 

At this stage, focus will move to optimisation of the code performance.

### Scope of the validation

This test validates the Tetracorder numerical evaluation and decision pipeline and reference-spectrum preparation against the native Tetracorder 6.00 Cuprite95 run.

The purpose of this test is to determine whether `python-tetracorder` reproduces the material classifications of native Tetracorder 6.00.

## Next steps

 - performance optimisation;  
   - refactoring  
   - numba jit where possible  
   - caching  

## Scope

This repository is **not** currently intended to be:

- a drop-in replacement for the complete Tetracorder software package;
- a reproduction of its graphical tools;
- a validated operational mineral-mapping product;
- an authoritative implementation of every historical Tetracorder behaviour;
- a substitute for expert spectral interpretation.

It is an independent Python reimplementation of the parts of the algorithm required to understand, test, and reproduce the published rule-based spectral evaluation workflow.

## Disclaimer

This software is provided for research and development purposes.

There is **no guarantee or warranty** that the implementation is:  
- scientifically equivalent to the original Tetracorder implementation, 
- complete, 
- correct,
- or suitable for any particular application.
  
Do not use results from this repository as the sole basis for scientific, commercial, exploration, environmental, safety, or other consequential decisions.

