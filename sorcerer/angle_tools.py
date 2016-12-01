import numpy as np


def bear(ra1, dec1, ra2, dec2):
    """
    Calculate the bearing of point b from point a.
    bearing is East of North [0,360)
    position angle is East of North (-180,180]
    """
    dlon = ra2 - ra1
    # dlat = dec2 - dec1
    y = np.sin(np.radians(dlon)) * np.cos(np.radians(dec2))
    x = np.cos(np.radians(dec1)) * np.sin(np.radians(dec2))
    x -= np.sin(np.radians(dec1)) * np.cos(np.radians(dec2)) * np.cos(np.radians(dlon))
    return np.degrees(np.arctan2(y, x))


# The following functions are explained
# at http://www.movable-type.co.uk/scripts/latlong.html
# phi ~ lat ~ Dec
# lambda ~ lon ~ RA
def gcd(ra1, dec1, ra2, dec2):
    """
    Great circle distance as calculated by the haversine formula
    ra/dec in degrees
    returns:
    sep in degrees
    """
    # TODO:  Vincenty formula
    # see - https://en.wikipedia.org/wiki/Great-circle_distance
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(np.radians(dlat) / 2) ** 2
    a += np.cos(np.radians(dec1)) * np.cos(np.radians(dec2)) * np.sin(np.radians(dlon) / 2) ** 2
    sep = np.degrees(2 * np.arcsin(min(1, np.sqrt(a))))
    return sep


def translate(ra, dec, r, theta):
    """
    Translate the point (ra,dec) a distance r (degrees)
    along angle theta (degrees)
    The translation is taken along an arc of a great circle.
    Return the (ra,dec) of the translated point.
    """
    factor = np.sin(np.radians(dec)) * np.cos(np.radians(r))
    factor += np.cos(np.radians(dec)) * np.sin(np.radians(r)) * np.cos(np.radians(theta))
    dec_out = np.degrees(np.arcsin(factor))

    y = np.sin(np.radians(theta)) * np.sin(np.radians(r)) * np.cos(np.radians(dec))
    x = np.cos(np.radians(r)) - np.sin(np.radians(dec)) * np.sin(np.radians(dec_out))
    ra_out = ra + np.degrees(np.arctan2(y, x))
    return ra_out, dec_out
