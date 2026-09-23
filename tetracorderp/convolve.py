"""
Vendored from emit-sds/tetracorder-lite (tetrapy package),
file tetracorder/conv/convolve.py 

Copyright (c) 2026 California Institute of Technology ("Caltech").
U.S. Government sponsorship acknowledged.
Licensed under the Apache License, Version 2.0.
See licenses/LICENSE-APACHE-2.0 and NOTICE.

MODIFIED FROM THE ORIGINAL: handling of specpr deleted points has been
removed; reference data is sanitised before it reaches this module.
########################################################################

Gaussian high-to-low resolution convolution (specpr function 17).

Tetracorder matches observed reflectance against a spectral library convolved to the
instrument's channels. Each native-resolution master spectrum is resampled onto the
target grid: for every output channel an effective Gaussian width is taken in
quadrature against the instrument's own resolution
(``sqrt(out_fwhm**2 - in_fwhm**2)``), a Gaussian bandpass is built on the native
grid, and the output sample is the normalized, channel-width-weighted average of the
input reflectance under that bandpass.

A :class:`Convolver` is constructed once for a target grid and reused for every
spectrum. The specpr "deleted point" sentinel is :data:`~tetrapy.conv.specpr.DELETED`;
values below :data:`DEL_THRESH` are treated as deleted.
"""

import numpy as np

#from tetrapy.conv.specpr import DELETED

DEL_THRESH = -1.0e30      # values below this are the deleted sentinel ## NOT USED HERE
FOUR_LN2 = 2.772589       # 4*ln(2); Gaussians are parameterised by FWHM
TLIM = 0.1e-7             # bandpass weights below this are ignored (f17 default)
EXP_FLOOR = -82.0         # exp() underflows to 0 below this exponent
MIN_WIDTH = 1.0e-10       # floor on the effective width to avoid divide-by-zero


class Convolver:
    """
    Resample native-resolution spectra onto a fixed target wavelength grid.

    This class implements specpr function 17 for Gaussian high-to-low resolution
    convolution. It constructs effective Gaussian bandpasses for each output channel
    by taking the instrument resolution in quadrature with the native resolution,
    then computes a weighted average of the input spectrum under each bandpass.

    Attributes
    ----------
    out_wave : numpy.ndarray
        Target wavelength grid (float64).
    out_fwhm : numpy.ndarray
        Target FWHM resolution at each wavelength (float64).
    """

    def __init__(self, out_wave: np.ndarray, out_fwhm: np.ndarray) -> None:
        """
        Initialize the convolver with a target wavelength/FWHM grid.

        Parameters
        ----------
        out_wave : numpy.ndarray
            Target output wavelengths (will be cast to float64).
        out_fwhm : numpy.ndarray
            Target output FWHM resolution at each wavelength (will be cast to float64).
            Must have the same length as ``out_wave``.
        """
        self.out_wave = np.asarray(out_wave, dtype=np.float64)
        self.out_fwhm = np.asarray(out_fwhm, dtype=np.float64)

    def convolve(self, in_wave: np.ndarray, in_fwhm: np.ndarray,
                 in_spec: np.ndarray) -> np.ndarray:
        """
        Convolve one native spectrum onto the target grid.

        Parameters
        ----------
        in_wave, in_fwhm, in_spec : numpy.ndarray
            The spectrum's native wavelengths, bandpass FWHM, and reflectance, all
            on the same native grid and the same length.

        Returns
        -------
        numpy.ndarray
            Convolved reflectance on the target grid (float32). Output channels with
            no valid support carry the deleted sentinel.
        """
        in_wave = np.asarray(in_wave, dtype=np.float64)
        in_fwhm = np.asarray(in_fwhm, dtype=np.float64)
        in_spec = np.asarray(in_spec, dtype=np.float64)
        if not (in_fwhm.size == in_wave.size == in_spec.size):
            raise ValueError("input wavelength/fwhm/spectrum length mismatch")

        # The deleted sentinel is no longer used as the sanitised spectra are passed 
        # from the SQLiteDB
        #valid = ~(_deleted(in_wave) | _deleted(in_spec))
        valid = np.ones(in_wave.size, dtype=bool)

        spacing = _channel_spacing(in_wave, valid)
        centers, widths, out_deleted = self._effective_gaussians(in_wave, in_fwhm)
        gauss = _bandpass(centers, widths, in_wave)

        # convol.r: drop weights below tlim or on deleted/invalid input channels
        weights = np.where(gauss < TLIM, 0.0, gauss)
        weights[:, ~valid] = 0.0

        contrib = weights * spacing[None, :]     # bandpass * channel width
        total = contrib.sum(axis=1)
        summed = contrib @ in_spec

        # The deleted sentinel is no loneger used, unsupported channels are NaN
        #out = np.full(self.out_wave.size, DELETED, dtype=np.float64)
        out = np.full(self.out_wave.size, np.nan, dtype=np.float64)
        good = (~out_deleted) & (total > 0.0) & (summed != 0.0)
        out[good] = summed[good] / total[good]
        return out.astype(np.float32)

    def _effective_gaussians(
        self, in_wave: np.ndarray, in_fwhm: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Per-output-channel Gaussian center, width, and deletion mask (gfiles.r).

        For each output channel the nearest input wavelength sets the native
        resolution subtracted in quadrature. If the output is no wider than the
        input the channel is degenerate: its width collapses to 1% of the native
        FWHM and its center snaps to the nearest input wavelength. Channels outside
        the native wavelength range are marked deleted.

        Parameters
        ----------
        in_wave : numpy.ndarray
            Native wavelength grid.
        in_fwhm : numpy.ndarray
            Native FWHM resolution at each wavelength.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
            ``(centers, widths, deleted)`` where ``centers`` are the effective
            Gaussian centers, ``widths`` are the effective widths (after quadrature
            subtraction), and ``deleted`` is a boolean mask of output channels
            outside the native wavelength range.
        """
        wav_min, wav_max = float(in_wave.min()), float(in_wave.max())
        centers = self.out_wave.copy()
        widths = self.out_fwhm.copy()
        deleted = np.zeros(self.out_wave.size, dtype=bool)

        for i, center in enumerate(centers):
            if not (wav_min < center < wav_max):
                deleted[i] = True
                continue
            near = int(np.argmin(np.abs(in_wave - center)))
            native = in_fwhm[near]
            if widths[i] <= native:
                widths[i] = native * 0.01
                centers[i] = in_wave[near]
            else:
                widths[i] = np.sqrt(widths[i] ** 2 - native ** 2)

        return centers, widths, deleted


# =======DELETED sentinel no longer used========================================
# def _deleted(a: np.ndarray) -> np.ndarray:
#     """
#     Check if array values represent deleted specpr data points.
# 
#     Parameters
#     ----------
#     a : numpy.ndarray
#         Array to check for deleted sentinel values.
# 
#     Returns
#     -------
#     numpy.ndarray
#         Boolean mask where True indicates deleted points (values < DEL_THRESH).
#     """
#     return a < DEL_THRESH
# =============================================================================


def _channel_spacing(wave: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """
    Local channel spacing for every input channel (delx.r).

    Interior channels use half the gap between their nearest valid neighbours on
    each side; edge channels use the one-sided gap. Channels with no valid
    neighbour get 0. This depends only on the native grid, so it is computed once
    per spectrum.

    Parameters
    ----------
    wave : numpy.ndarray
        Wavelength values for each channel.
    valid : numpy.ndarray
        Boolean mask indicating which channels are valid (not deleted).

    Returns
    -------
    numpy.ndarray
        Local spacing for each channel, used to weight contributions in the
        convolution. Shape matches ``wave``.
    """
    n = wave.size
    spacing = np.zeros(n, dtype=np.float64)
    idx = np.flatnonzero(valid)
    if n < 2 or idx.size == 0:
        return spacing

    for i in range(n):
        above = idx[np.searchsorted(idx, i, side="right"):]
        below = idx[:np.searchsorted(idx, i, side="left")]
        hi = above[0] if above.size else -1
        lo = below[-1] if below.size else -1

        if hi >= 0 and lo >= 0:
            spacing[i] = abs(wave[hi] - wave[lo]) * 0.5
        elif hi >= 0:
            spacing[i] = abs(wave[i] - wave[hi])
        elif lo >= 0:
            spacing[i] = abs(wave[i] - wave[lo])

    return spacing


def _bandpass(
    centers: np.ndarray, widths: np.ndarray, in_wave: np.ndarray
) -> np.ndarray:
    """
    Gaussian bandpass matrix ``(n_out, n_in)`` (ggauss.r).

    Computes a Gaussian weight matrix where each row is the bandpass for one
    output channel. The Gaussian form is
    ``exp(-4 ln2 / width**2 * (in_wave - center)**2)``, with deleted input
    wavelengths and underflowing exponents set to zero. If every weight for an
    output channel underflows, its single nearest non-deleted input channel is
    given weight 1.0 -- how the shipped library extrapolates channels beyond a
    spectrum's native coverage.

    Parameters
    ----------
    centers : numpy.ndarray
        Gaussian center wavelength for each output channel, shape ``(n_out,)``.
    widths : numpy.ndarray
        Gaussian FWHM width for each output channel, shape ``(n_out,)``.
    in_wave : numpy.ndarray
        Input wavelength grid, shape ``(n_in,)``.

    Returns
    -------
    numpy.ndarray
        Gaussian weight matrix, shape ``(n_out, n_in)``, where ``[i, j]`` is
        the contribution weight of input channel j to output channel i.
    """
    # remove the "deleted" infrastructure
    #del_in = _deleted(in_wave)
    del_in = np.zeros(in_wave.size, dtype=bool)
    safe_widths = np.where(widths <= MIN_WIDTH, MIN_WIDTH, widths)
    diff = in_wave[None, :] - centers[:, None]
    exponent = (-FOUR_LN2 / safe_widths[:, None] ** 2) * diff ** 2

    gauss = np.where(exponent < EXP_FLOOR, 0.0, np.exp(exponent))
    gauss[:, del_in] = 0.0

    all_zero = ~gauss.any(axis=1)
    if all_zero.any():
        nearest = np.where(del_in[None, :], np.inf, np.abs(diff)).argmin(axis=1)
        for i in np.flatnonzero(all_zero):
            gauss[i, nearest[i]] = 1.0

    return gauss
