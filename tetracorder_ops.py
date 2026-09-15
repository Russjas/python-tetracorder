"""tetracorder operations reimplemented in numpy scipy.
to be used by the rule interpreter algorithm once implemented"""

from dataclasses import dataclass
from typing import Literal
from copy import deepcopy
import re

import numpy as np
from scipy.interpolate import CubicSpline



#=============== pre-defined values for the symbolic thresholds in Tetracorder
#comments are preserved from Tetracorder cmd files
##comments are mine
DEFAULT_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
#GLBLFITVEGR is for vegetation red-edge detection, case 1
"[GLBLFITVEGR]": [0.3, 0.40],
#GLBLFITVEGT is for vegetation type detection, case 2
"[GLBLFITVEGT]": [0.2, 0.30],
#GLBLFITVEGT is for vegetation water detection, case 3, 4, 5
"[GLBLFITVEGW]": [0.4, 0.55],
#GLBLDPFITREE is set high to reduce false positives
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

#  Continuum thresholds

#NB: ct, rct, lct do not use fuzzy logic
"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

##fuzzy logic params for feature tests

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
}
DEFAULT_BAK23_VALUES = {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
"[GLBLFITVEGR]": [0.3, 0.40],
"[GLBLFITVEGT]": [0.2, 0.30],
"[GLBLFITVEGW]": [0.4, 0.55],
"[GLBLDPFITREE]": [0.65, 0.80],## These are ommited from the bak23 variable cmd file - added here to avoid leaving symbols in the rule
"[GLBLFITREE]": [0.5, 0.6],## These are ommited from the bak23 variable cmd file - added here to avoid leaving symbols in the rule

"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.001, 0.002],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
}
EMIT_C_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
"[GLBLFITVEGR]": [0.3, 0.40],
"[GLBLFITVEGT]": [0.2, 0.30],
"[GLBLFITVEGW]": [0.4, 0.55],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
}
MMM_09C_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.005,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.01,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.02,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.04,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
}
MMM255T_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.005,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.01,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.02,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.04,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
}
HYB2RYUG_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.001,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.002,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.002,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.0025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.004,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.0002,  0.0004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.00002, 0.00004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.0002, 0.0004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.0001, 0.00002],     # r*bd for 1-      micron region ##sic Tetracorder but right<left could be typo
"[RBD14]": [0.00008, 0.0002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.00007, 0.00011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.00007, 0.00011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.0001, 0.0002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.00007, 0.00011],    # r*bd for 2.5-    micron region #sic Tetracorder
"[RBD30]": [0.00005, 0.0001],     # r*bd for 3-      micron region
"[RBD35]": [0.00005, 0.0001],     # r*bd for 3.5-    micron region
}

#============== Preparing the rulesets for use by RuleEvaluator================

VARIABLE_PRESETS = {
    "default": DEFAULT_VALUES,
    "default_bak23": DEFAULT_BAK23_VALUES,
    "emit_c": EMIT_C_VALUES,
    "MMM_09c": MMM_09C_VALUES,
    "MMM255t": MMM255T_VALUES,
    "HYB2ryug": HYB2RYUG_VALUES,
}

CONSTRAINT_RE = re.compile(
    r"(FITALL>|FIT>|DEPTHALL>|DEPTH-FIT>|FD-FIT>|FDALL>|FD-DEPTH>|DEPTH>|FD>)"
    r"\s*"
    r"(\[[^\]]+\]|[-+]?\d*\.?\d+(?:\s+[-+]?\d*\.?\d+)?)"
)


def parse_constraints(constraints):
    parsed = []

    for line in constraints:
        line = line.removeprefix("constraint:").strip()

        for test, values in CONSTRAINT_RE.findall(line):
            parsed.append({
                "test": test,
                "values": values.strip(),
            })

    return parsed

def prepare_rules(rules, mode="default"):
    """
    Expects a rule from the parsed json rules json
    mode = default, default_bak23, emit_c, MMM_09c, MMM255t,
    HYB2ryug defined in Tetracorder 6.00a VARIABLES/cmd.lib.setup.variables-*
    Custom modes can be established by creating a new preset dictionary
    """
    values = VARIABLE_PRESETS[mode]

    
    def replace(value):

        if isinstance(value, dict):
            return {key: replace(item) for key, item in value.items()}

        if isinstance(value, list):
            # Parser stored symbolic values as one-item lists.
            # Collapse the list and substitute the symbol directly.
            if (len(value) == 1 and isinstance(value[0], str) and value[0] in values):
                return deepcopy(values[value[0]])

            return [replace(item) for item in value]

        if isinstance(value, str) and value in values:
            return deepcopy(values[value])
        

        return value

    return replace(rules)

#==========================================================================

class GaussianConvolver:
    """
    Convolve high-resolution library spectra to scanner channels using
    Gaussian spectral response functions defined by scanner FWHM.

    The convolution weights are calculated once and can then be reused
    for multiple library spectra on the same wavelength grid.

    Usage: 
    convolver = GaussianConvolver(
                lib_wl=library_wavelengths,
                scanner_wl=scanner_wavelengths,
                scanner_fwhm=scanner_fwhm,
                )
    convolved_reflectance = convolver.convolve(library_reflectance,)
    """

    def __init__(self, lib_wl, scanner_wl, scanner_fwhm):
        self.lib_wl = np.asarray(lib_wl, dtype=float)
        self.scanner_wl = np.asarray(scanner_wl, dtype=float)
        self.scanner_fwhm = np.asarray(scanner_fwhm, dtype=float)

        if self.scanner_fwhm.ndim == 0 or self.scanner_fwhm.size == 1:
            self.scanner_fwhm = np.full(
                self.scanner_wl.shape,
                float(self.scanner_fwhm.ravel()[0]),
            )

        if self.scanner_wl.shape != self.scanner_fwhm.shape:
            raise ValueError(
                "scanner_wl and scanner_fwhm must have equal lengths"
            )

        if self.lib_wl.ndim != 1:
            raise ValueError("lib_wl must be one-dimensional")

        if np.any(np.diff(self.lib_wl) <= 0):
            raise ValueError("lib_wl must be strictly increasing")

        if np.any(self.scanner_fwhm <= 0):
            raise ValueError("All scanner FWHM values must be positive")

        self.weights = self._build_weights()

    def _build_weights(self):
        # Gaussian response:
        # exp(-4 ln(2) * ((x - centre) / FWHM)**2)
        offset = (
            self.lib_wl[None, :]
            - self.scanner_wl[:, None]
        )

        weights = np.exp(
            -4.0 * np.log(2.0)
            * (offset / self.scanner_fwhm[:, None]) ** 2
        )

        # Approximate integration widths on the library wavelength grid.
        # This matters if the library grid is not perfectly uniform.
        spacing = np.gradient(self.lib_wl)
        weights *= spacing[None, :]

        # Do not calculate scanner bands whose centres fall outside the
        # library wavelength coverage.
        outside = (
            (self.scanner_wl <= self.lib_wl.min())
            | (self.scanner_wl >= self.lib_wl.max())
        )
        weights[outside] = 0.0

        # Negligible weights add computation without affecting the result.
        weights[weights < 1e-8] = 0.0

        return weights

    def convolve(self, lib_refl):
        """
        Parameters
        ----------
        lib_refl : ndarray
            Either:
                (library_bands,)
            or:
                (number_of_spectra, library_bands)

        Returns
        -------
        ndarray
            Spectrum/spectra convolved to scanner_wl.
        """
        spectra = np.asarray(lib_refl, dtype=float)
        was_1d = spectra.ndim == 1

        if was_1d:
            spectra = spectra[None, :]

        if spectra.ndim != 2:
            raise ValueError("lib_refl must be one- or two-dimensional")

        if spectra.shape[1] != self.lib_wl.size:
            raise ValueError(
                f"Expected {self.lib_wl.size} library bands, "
                f"received {spectra.shape[1]}"
            )

        valid = np.isfinite(spectra)

        # Renormalise independently where individual spectra contain NaNs.
        numerator = (
            np.where(valid, spectra, 0.0)
            @ self.weights.T
        )

        denominator = (
            valid.astype(float)
            @ self.weights.T
        )

        result = np.full(numerator.shape, np.nan, dtype=float)

        np.divide(
            numerator,
            denominator,
            out=result,
            where=denominator > 0,
        )

        if was_1d:
            return result[0]

        return result


@dataclass
class ContinuumFeature:
    """A class for holding a continuum removed feature for the python tetracorder algorithm"""
    continuum_type: Literal["linear", "curved"]

    wavelengths: np.ndarray
    continuum_removed: np.ndarray
    continuum: np.ndarray

    left_continuum: np.ndarray
    right_continuum: np.ndarray
# Validation for wavelenth dataset

def _validate_wavelengths(wavelengths):
    """
    Validate and return a one-dimensional, increasing wavelength array.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)

    if wavelengths.ndim != 1:
        raise ValueError(
            f"Wavelengths must be one-dimensional; got shape "
            f"{wavelengths.shape}."
        )

    if wavelengths.size == 0:
        raise ValueError("Wavelength array is empty.")

    if not np.all(np.isfinite(wavelengths)):
        raise ValueError("Wavelength array contains non-finite values.")

    if not np.all(np.diff(wavelengths) > 0):
        raise ValueError(
            "Wavelengths must be strictly increasing with no duplicates."
        )

    return wavelengths


def _validate_window(wavelengths, window, name="window"):
    """
    Check that a window contains two valid bounds and at least one band.

    Returns
    -------
    valid : bool
    message : str
    """
    try:
        lower, upper = window
    except (TypeError, ValueError):
        return (
            False,
            f"{name} must contain exactly two wavelength bounds; "
            f"got {window!r}.",
        )

    try:
        lower = float(lower)
        upper = float(upper)
    except (TypeError, ValueError):
        return (
            False,
            f"{name} bounds must be numeric; got {window!r}.",
        )

    if not np.isfinite(lower) or not np.isfinite(upper):
        return (
            False,
            f"{name} contains a non-finite wavelength bound: "
            f"{window!r}.",
        )

    if lower >= upper:
        return (
            False,
            f"{name} lower bound must be less than its upper bound; "
            f"got {lower:g}–{upper:g} nm.",
        )

    mask = (wavelengths >= lower) & (wavelengths <= upper)

    if not np.any(mask):
        return (
            False,
            f"{name} {lower:g}–{upper:g} nm contains no bands.",
        )

    return True, ""



# Continuum calculation ported from specspr


def linear_feature_continuum(
    spectra,
    wavelengths,
    left_window,
    right_window,
):
    """
    Reproduce the linear continuum used by Specpr/Tetracorder.

    spectra:
        (..., bands)

    wavelengths:
        (bands,)

    Returns
    -------
    continuum:
        (..., feature_bands)

    continuum_removed:
        (..., feature_bands)

    feature_wavelengths:
        (feature_bands,)

    left_mean, right_mean:
        Mean continuum reflectances.
    """
    spectra = np.asanyarray(spectra)
    wavelengths = _validate_wavelengths(wavelengths)

    if spectra.shape[-1] != wavelengths.size:
        raise ValueError(
            f"Spectra have {spectra.shape[-1]} bands but "
            f"{wavelengths.size} wavelengths were supplied."
        )

    check, msg = _validate_window(
        wavelengths,
        left_window,
        name="Left continuum window",
    )
    if not check:
        raise ValueError(msg)

    check, msg = _validate_window(
        wavelengths,
        right_window,
        name="Right continuum window",
    )
    if not check:
        raise ValueError(msg)
    left = (
        (wavelengths >= left_window[0])
        & (wavelengths <= left_window[1])
    )

    right = (
        (wavelengths >= right_window[0])
        & (wavelengths <= right_window[1])
    )

    feature = (
        (wavelengths > left_window[1])
        & (wavelengths < right_window[0])
    )

    if not feature.any():
        raise ValueError("No bands in feature interval")

    left_refl = spectra[..., left].mean(axis=-1)
    right_refl = spectra[..., right].mean(axis=-1)

    # Specpr averages the actual channel wavelengths.
    left_wl = wavelengths[left].mean()
    right_wl = wavelengths[right].mean()

    slope = (
        (right_refl - left_refl)
        / (right_wl - left_wl)
    )

    intercept = right_refl - slope * right_wl

    feature_wl = wavelengths[feature]

    continuum = (
        slope[..., None] * feature_wl
        + intercept[..., None]
    )

    continuum_removed = spectra[..., feature] / continuum

    return ContinuumFeature(
    continuum_type="linear",
    wavelengths=feature_wl,
    continuum_removed=continuum_removed,
    continuum=continuum,
    left_continuum=left_refl,
    right_continuum=right_refl,
)



def curved_feature_continuum(
    spectra,
    wavelengths,
    continuum_windows,
):
    """
    Reproduce the curved continuum used by Specpr/Tetracorder.

    continuum_windows:
        [
            (outer_left_start, outer_left_stop),
            (inner_left_start, inner_left_stop),
            (inner_right_start, inner_right_stop),
            (outer_right_start, outer_right_stop),
        ]
    """
    spectra = np.asanyarray(spectra)
    wavelengths = _validate_wavelengths(wavelengths)

    if spectra.shape[-1] != wavelengths.size:
        raise ValueError(
            f"Spectra have {spectra.shape[-1]} bands but "
            f"{wavelengths.size} wavelengths were supplied."
        )

    if len(continuum_windows) != 4:
        raise ValueError(
            "Curved continuum requires four wavelength windows"
        )

    for i, window in enumerate(continuum_windows, start=1):
        check, msg = _validate_window(
            wavelengths,
            window,
            name=f"Continuum window {i}",
        )
        if not check:
            raise ValueError(msg)


    anchor_wavelengths = np.array(
        [
            (start + stop) / 2.0
            for start, stop in continuum_windows
        ],
        dtype=float,
    )

    anchor_reflectances = []

    for start, stop in continuum_windows:
        selected = (
            (wavelengths >= start)
            & (wavelengths <= stop)
        )

        if not selected.any():
            raise ValueError(
                f"No bands in continuum window {start}–{stop}"
            )

        anchor_reflectances.append(
            spectra[..., selected].mean(axis=-1)
        )

    # Shape becomes (..., 4)
    anchor_reflectances = np.stack(
        anchor_reflectances,
        axis=-1,
    )

    feature = (
        (wavelengths > continuum_windows[1][1])
        & (wavelengths < continuum_windows[2][0])
    )

    if not feature.any():
        raise ValueError("No bands in feature interval")

    feature_wl = wavelengths[feature]

    spline = CubicSpline(
        anchor_wavelengths,
        anchor_reflectances,
        axis=-1,
        bc_type="natural",
    )

    continuum = spline(feature_wl)
    continuum_removed = spectra[..., feature] / continuum

    return ContinuumFeature(
    continuum_type="curved",
    wavelengths=feature_wl,
    continuum_removed=continuum_removed,
    continuum=continuum,
    left_continuum=anchor_reflectances[..., 1],
    right_continuum=anchor_reflectances[..., 2],
)

# ===============Characterise primitive ========================
def characterise_feature(reference, target, polarity = "absorption"):
    """
    Characterise the relationship between a reference feature and one
    or more target features using the Specpr/Tetracorder fitting method.

    Parameters
    ----------
    reference : ContinuumFeature
        A single reference feature. Its continuum_removed array must
        have shape (bands,).

    target : ContinuumFeature
        Target feature or features. Its continuum_removed array may
        have shape (bands,) or (..., bands).

    Returns
    -------
    dict
        Feature fit, depth and continuum-shape measurements.
    """
    if polarity not in {"absorption", "emission"}:
        raise ValueError("Polarity must be 'absorption' or 'emission'")

    if reference.continuum_type != target.continuum_type:
        raise ValueError("Reference and target use different continuum types")

    reference_wl = np.asarray(reference.wavelengths, dtype=float)
    target_wl = np.asarray(target.wavelengths, dtype=float)

    if (reference_wl.shape != target_wl.shape
        or not np.allclose(reference_wl, target_wl)):
        raise ValueError("Reference and target feature wavelengths do not match.")

    reference_cr = np.ma.asarray(reference.continuum_removed, dtype=float,)

    target_cr = np.ma.asarray(target.continuum_removed, dtype=float,)

    if reference_cr.ndim != 1:
        raise ValueError("Reference feature must be one-dimensional.")

    if target_cr.shape[-1] != reference_cr.size:
        raise ValueError("Reference and target features contain different numbers of bands.")
    #2 is set as a mathematically acurate threshold, 
    # but may increase to screen ridiculously narrow features for the sensor
    if reference_cr.count() < 2: 
        raise ValueError("Reference feature has fewer than two valid bands.")

    # --------------------------------------------------------------
    # Determine whether the exemplar is primarily an absorption or
    # emission feature.
    # --------------------------------------------------------------
    if polarity == "absorption":
        extremum_index = int(np.ma.argmin(reference_cr))
    elif polarity == "emission":
        extremum_index = int(np.ma.argmax(reference_cr))

    reference_depth = (1.0 - reference_cr[extremum_index])

    # reference_depth is positive for absorption and negative for
    # emission.

    # --------------------------------------------------------------
    # Mean-centred reference and target feature shapes.
    # --------------------------------------------------------------

    reference_deviation = (reference_cr - reference_cr.mean())

    target_deviation = (target_cr - target_cr.mean(axis=-1, keepdims=True))

    # tetracorder uses the name covariance, so preserved here
    covariance = np.ma.sum(target_deviation * reference_deviation, axis=-1,)

    reference_variance = np.ma.sum(reference_deviation**2)

    target_variance = np.ma.sum(target_deviation**2, axis=-1,)

    if (np.ma.is_masked(reference_variance) or reference_variance <= 0):
        raise ValueError("Reference feature has no measurable variation.")

    # --------------------------------------------------------------
    # Scale reference-feature depth to the target.
    # --------------------------------------------------------------

    slope = covariance / reference_variance

    # --------------------------------------------------------------
    # Specpr goodness of fit.
    #
    # This is equivalent to the absolute Pearson correlation, but is
    # retained in the form used by the original fitting calculation.
    # --------------------------------------------------------------

    with np.errstate(divide="ignore", invalid="ignore"):
        reverse_slope = covariance / target_variance
        fit = np.ma.sqrt(np.ma.abs(slope * reverse_slope))

    invalid_fit = (
            (np.ma.filled(reference_variance, 0.0) <= 0)
            | (np.ma.filled(target_variance, 0.0) <= 0)
            | (np.abs(np.ma.filled(slope, 0.0)) < 1e-21)
        )
    
    fit = np.ma.masked_where(invalid_fit, fit)
    fit = np.ma.clip(fit, 0.0, 1.0)

    # No need to construct the full fitted-reference cube:
    # fitted_reference = 1 + slope * (reference_cr - 1)
    # Its depth at the reference extremum simplifies to this.

    depth = slope * reference_depth
    fit_depth = fit * depth

    # --------------------------------------------------------------
    # Target continuum geometry at the reference extremum.
    # --------------------------------------------------------------

    continuum = np.ma.asarray(target.continuum, dtype=float,)[..., extremum_index]

    left_continuum = np.ma.asarray(target.left_continuum, dtype=float,)

    right_continuum = np.ma.asarray(target.right_continuum, dtype=float,)

    band_bottom = continuum * (1.0 - depth)

    # tetracorder returns absolute values here after infering polarity
    # in this implementation polarity is provided thus signed value
    # is returned.
    reflectance_depth = continuum * depth

    # --------------------------------------------------------------
    # Continuum and shoulder ratios.
    # --------------------------------------------------------------
    # Tetracorder implements a 1e-7 safe denominator. I do not implement this here
    #  Non-finite results are retained so the rule evaluator can reject them.
    with np.errstate(divide="ignore", invalid="ignore"):
        left_right_ratio = (left_continuum / right_continuum)

        right_left_ratio = (right_continuum / left_continuum)

        left_shoulder_ratio = ((left_continuum - band_bottom) / (right_continuum - band_bottom))

        right_shoulder_ratio = ((right_continuum - band_bottom) / (left_continuum - band_bottom))

    return {
        "fit": fit,
        "depth": depth,
        "fit_depth": fit_depth,
        "slope": slope,

        "polarity": polarity,
        "extremum_index": extremum_index,
        "reference_depth": reference_depth,

        "continuum": continuum,
        "left_continuum": left_continuum,
        "right_continuum": right_continuum,
        "band_bottom": band_bottom,
        "reflectance_depth": reflectance_depth,

        "left_right_ratio": left_right_ratio,
        "right_left_ratio": right_left_ratio,
        "left_shoulder_ratio": left_shoulder_ratio,
        "right_shoulder_ratio": right_shoulder_ratio,

        "value_count": reference_cr.count(),
        "polarity": polarity,
    }


#==================== 

