import math
import numpy as np


def draw_ellipse(X, Y, x, y, major, minor, pa):
    """
    Returns boolean array if (x,y) from X,Y is located in or on the ellipse

    Args:
        X, Y: numpy arrays of equal length, from which we generate
              cartesian coordinate pairs
        major, minor: The elliptical major and minor axes, in the
                      units as the X,Y points
        pa: The angle (degrees) the major axis makes to the x-axis,
            going counter clockwise.
    """
    # Perform translation and rotation first
    Xdash = (X-x) * np.cos(np.radians(pa)) + (Y-y) * np.sin(np.radians(pa))
    Ydash = (Y-y) * np.cos(np.radians(pa)) - (X-x) * np.sin(np.radians(pa))
    return (Xdash**2 / major**2 + (Ydash)**2 / minor**2 <= 1)


def elliptical_gaussian(X, Y, x, y, amp, major, minor, pa):
    """
    Generate a model 2d Gaussian with the given parameters.
    Evaluate this model at the given locations x,y.

    :param x,y: locations at which to calculate values
    :param xo,yo: position of Gaussian
    :param amp: amplitude of Gaussian
    :param major,minor: axes (sigmas)
    :param theta: position angle (degrees) CCW from x-axis
    :return: Gaussian function evaluated at x,y locations
    """
    try:
        sint, cost = math.sin(np.radians(pa)), math.cos(np.radians(pa))
    except ValueError as e:
        if 'math domain error' in e.args:
            sint, cost = np.nan, np.nan
        else:
            raise
    xxo = X - x
    yyo = Y - y
    exp = (xxo * cost + yyo * sint) ** 2 / major ** 2 + \
          (xxo * sint - yyo * cost) ** 2 / minor ** 2
    exp *= -1. / 2
    return amp * np.exp(exp)


def draw_gaussian(source, data):
    # Only apply Gaussian's out to this limit. Dramatically speeds up image creation.
    gaussian_limit = 6 * max(source.major, source.minor)

    ysize, xsize = data.shape
    xmin = int(max(0, source.loc_x - gaussian_limit))
    xmax = int(min(xsize, source.loc_x + gaussian_limit))
    ymin = int(max(0, source.loc_y - gaussian_limit))
    ymax = int(min(ysize, source.loc_y + gaussian_limit))
    x, y = np.mgrid[xmin:xmax, ymin:ymax]
    x = x.ravel()
    y = y.ravel()

    model = elliptical_gaussian(x, y, source.loc_x, source.loc_y,
                                source.peak, source.major,
                                source.minor, source.pa)
    data[y, x] += model