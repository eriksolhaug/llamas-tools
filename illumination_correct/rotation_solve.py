'''
Rotation solve utilities.
'''

import numpy as np
from scipy.ndimage import rotate, median_filter
from scipy.optimize import minimize_scalar

def whitelight(cube, wmin=None, wmax=None):
    """
    Generate a white light image using the median of the data cube.

    cube : (nwave, ny, nx) array, NaN where there is no data.
    wmin, wmax : plane INDICES (not Angstrom) of the wavelength range. Use
        np.searchsorted(wave, wavelength) to convert. None = full range.
    """
    if wmin is not None and wmax is not None:
        cube = cube[wmin:wmax, :, :]

    # Compute median along the wavelength axis, ignoring the NaNs of the empty spaxels
    white_light_image = np.nanmedian(cube, axis=0)

    return white_light_image

def make_footprint(shape, angle):
    """
    Create a footprint for the median filter based on the specified shape and rotation angle.
    """
    nx, ny = shape

    # Odd square canvas, large enough to hold the box at any angle (its diagonal)
    side = int(np.ceil(np.hypot(nx, ny)))
    side += 1 - side % 2
    center = side // 2

    # Fill the box in the middle of the canvas with 1s
    canvas = np.zeros((side, side))
    y0 = center - ny // 2
    x0 = center - nx // 2
    canvas[y0:y0 + ny, x0:x0 + nx] = 1

    # Rotate; order=0 keeps the values 0 or 1. scipy rotates clockwise for the image
    # displayed with origin='lower', so the angle is negated to make it counterclockwise
    footprint = rotate(canvas, -angle, reshape=False, order=0)

    return footprint.astype(bool)

def residual_rms(image, shape, angle):
    """
    RMS of the image minus its median filtered version, for a footprint of the given
    shape and angle. It is smallest when the footprint lies along the stripes, because
    the filter then keeps the stripe pattern instead of averaging across it.

    The NaNs are filled with the median of the image before filtering and are left out
    of the RMS.
    """
    valid = np.isfinite(image)
    filled = np.where(valid, image, np.nanmedian(image))
    footprint = make_footprint(shape, angle)
    filtered = median_filter(filled, footprint=footprint, mode="nearest")
    return np.std((image - filtered)[valid])

def solve_rotation(image, shape, step=1.0, tol=0.05):
    """
    Fit the rotation angle of the footprint that best aligns with the pattern of the image.

    image : 2D white light image (NaN where there is no data).
    shape : (nx, ny) of the footprint, see make_footprint.
    step : angle spacing in degrees of the coarse scan used to bracket the minimum.
    tol : precision in degrees of the final fit.

    The angle is continuous: a coarse scan over [0, 180) finds the neighbourhood of the
    minimum of residual_rms, then scipy's bounded scalar minimizer refines it. The box
    is symmetric under a 180 degree turn, so the angle is periodic.

    Returns the angle in [0, 180) and the contrast of the scan, (max - min) / max of the
    coarse residual RMS, which is small when the image has no clear pattern to align with.
    """
    angles = np.arange(0, 180, step)
    scores = np.array([residual_rms(image, shape, a) for a in angles])
    best = angles[scores.argmin()]
    contrast = (scores.max() - scores.min()) / scores.max()

    # Fit within one step of the coarse minimum; it may cross 0 or 180, handled by the modulo
    fit = minimize_scalar(lambda a: residual_rms(image, shape, a),
                          bounds=(best - step, best + step), method="bounded",
                          options={"xatol": tol})

    return fit.x % 180, contrast
