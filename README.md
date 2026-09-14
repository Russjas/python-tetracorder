# python-tetracorder
A port of the rules and evaluation logic of the tetracorder algorithm to pure python - with no guarantee or waranty!

## Status

**Early development / experimental.**

This repository is currently a working area for reimplementing parts of the Tetracorder spectral identification workflow in Python.

The immediate aim is not to reproduce the full original software environment, user interface, or processing system. The focus is on the underlying spectral rules and evaluation logic: taking spectral reference data, applying Tetracorder-style rules, and reproducing the decision process in a form that is readable, testable, and usable from Python.

The implementation is incomplete and should not currently be treated as a validated replacement for Tetracorder.

## Goals

The project is intended to:

- port Tetracorder mineral identification rules into a structured, machine-readable form;
- reproduce the relevant spectral evaluation logic in pure Python;
- The implementation is therefore intended to rely primarily on standard Python scientific tools such as:
  - NumPy
  - SciPy
  - SQLite
- make the rule evaluation process easier to inspect, test, modify, and integrate into other Python workflows;
- provide a foundation for testing Tetracorder-style mineral identification against modern hyperspectral datasets.

## Rules

The project is based on published/released [Tetracorder](https://github.com/PSI-edu/spectroscopy-tetracorder) rule definitions and attempts to reproduce their evaluation behaviour independently in Python.
The first step is to parse these rules from the original cmd files into a json schema

## Reference spectra

Tetracorder rules reference spectra from USGS spectral libraries. The most recent version of the USGS spectral library is available [here](). However the tetracorder rules refer to specific refereence spectra used in the tetracorder rulesets are from a previous version of this library, and mapping is not straightforward.
The sql reference spectra database in this repository is derived from the libraries in the [Tetracorder repo](https://github.com/PSI-edu/spectroscopy-tetracorder). The database here is derived from the native libraries rather than the convolved variants.

## Evaluation logic  

in progress

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

There is **no guarantee or warranty** that the implementation is complete, correct, scientifically equivalent to the original Tetracorder implementation, or suitable for any particular application.

Do not use results from this repository as the sole basis for scientific, commercial, exploration, environmental, safety, or other consequential decisions.

