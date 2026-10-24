'''
Median filter for illumination correction.
'''

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")
from copy import deepcopy

def median_2dfilter(cube, footprint):
    """
    Median filter every wavelength plane of a cube with a 2D footprint, ignoring the NaNs.

    cube : (nwave, ny, nx) array, NaN where there is no data.
    footprint : square boolean array with an odd side, see rotation_solve.make_footprint.
        Its axes are (y, x), like the planes of the cube.

    The median of each pixel is taken over the valid pixels under the footprint only, so
    the empty edges of the field do not bias it (scipy's median_filter does not handle
    NaNs). The NaN pixels of the cube stay NaN in the output.

    Returns the median filtered cube, same shape as cube.
    """
    if footprint.shape[0] != footprint.shape[1] or footprint.shape[0] % 2 == 0:
        raise ValueError("footprint must be a square array with an odd side.")

    ny, nx = cube.shape[1:]
    pad = footprint.shape[0] // 2
    width = nx + 2 * pad

    # Offsets of the footprint pixels from its center, as steps in the flattened padded plane
    fy, fx = np.nonzero(footprint)
    offsets = (fy - pad) * width + (fx - pad)

    median_cube = np.full(cube.shape, np.nan, dtype=cube.dtype)
    padded = np.full((ny + 2 * pad, width), np.nan, dtype=cube.dtype)

    for k, plane in enumerate(cube):
        valid = np.isfinite(plane)
        if not valid.any():
            continue

        # NaN border so that the footprint can hang over the edge of the plane
        padded[pad:pad + ny, pad:pad + nx] = plane

        # For each valid pixel, the values under the footprint: (footprint pixels, valid pixels)
        yy, xx = np.nonzero(valid)
        centers = (yy + pad) * width + (xx + pad)
        values = padded.ravel()[centers[None, :] + offsets[:, None]]

        median_cube[k][valid] = np.nanmedian(values, axis=0)

    return median_cube


