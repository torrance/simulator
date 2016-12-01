from sourcerer.angle_tools import gcd, translate, bear
from astropy import wcs
import numpy as np


class WCSHelper:
    def __init__(self, header):
        self.wcs = wcs.WCS(header)
        self.header = header

    def pix2sky(self, pixel):
        """
        Take pixel=(x,y) coords
        convert to pos=(ra,dec) coords
        """
        x, y = pixel
        return self.wcs.wcs_pix2world([[x, y]], 1)[0]

    def sky2pix(self, pos):
        """
        Take pos = (ra,dec) coords
        convert to pixel = (x,y) coords
        """
        pixel = self.wcs.wcs_world2pix([pos], 1)
        return [pixel[0][0], pixel[0][1]]

    def pix2sky_ellipse(self, pixel, sx, sy, theta):
        """
        Convert an ellipse from pixel to sky coords
        sx/sy vectors are calculated at an origin pos=(x,y)
        Input parameters are:
        x,y - the x,y pixels corresponding to the ra/dec position
        sx, sy - the major minor axes (FWHM) in pixels
        theta - the position angle in degrees
        Output params are all in degrees

        :param pixel: [x,y] of the ellipse center
        :param sx: major axis
        :param sy: minor axis
        :param theta: position angle
        :return: ra, dec, a, b, pa
        """
        ra, dec = self.pix2sky(pixel)
        x, y = pixel
        v_sx = [x + sx * np.cos(np.radians(theta)),
                y + sx * np.sin(np.radians(theta))]
        ra2, dec2 = self.pix2sky(v_sx)
        major = gcd(ra, dec, ra2, dec2)
        pa = bear(ra, dec, ra2, dec2)

        v_sy = [x + sy * np.cos(np.radians(theta-90)),
                y + sy * np.sin(np.radians(theta-90))]
        ra2, dec2 = self.pix2sky(v_sy)
        minor = gcd(ra, dec, ra2, dec2)
        pa2 = bear(ra, dec, ra2, dec2) - 90

        # The a/b vectors are perpendicular in sky space,
        # but not always in pixel space
        # so we have to account for this by calculating
        # the angle between the two vectors
        # and modifying the minor axis length
        defect = pa - pa2
        minor *= abs(np.cos(np.radians(defect)))
        return ra, dec, major, minor, pa

    def sky2pix_ellipse(self, pos, a, b, pa):
        """
        Convert an ellipse from sky to pixel corrds
        a/b vectors are calculated at an origin pos=(ra,dec)
        All input parameters are in degrees
        Output parameters are:
        x,y - the x,y pixels corresponding to the ra/dec position
        sx, sy - the major minor axes (FWHM) in pixels
        theta - the position angle in degrees

        :param pos: [ra,dec] of the ellipse center
        :param a: major axis
        :param b: minor axis
        :param pa: position angle
        :return: x, y, sx, sy, theta
        """
        ra, dec = pos
        x, y = self.sky2pix(pos)

        x_off, y_off = self.sky2pix(translate(ra, dec, a, pa))
        sx = np.hypot((x - x_off), (y - y_off))
        theta = np.arctan2((y_off - y), (x_off - x))

        x_off, y_off = self.sky2pix(translate(ra, dec, b, pa-90))
        sy = np.hypot((x - x_off), (y - y_off))
        theta2 = np.arctan2((y_off - y), (x_off - x)) - np.pi/2

        # The a/b vectors are perpendicular in sky space,
        # but not always in pixel space
        # so we have to account for this by calculating
        # the angle between the two vectors
        # and modifying the minor axis length
        defect = theta - theta2
        sy *= abs(np.cos(defect))

        return x, y, sx, sy, np.degrees(theta)

    def beamarea_pix(self):
        """
        Returns the area of the beam in pixels squared,
        taking into account the special 4 ln(2) factor.
        """
        beamsigma1 = self.header['BMAJ'] / self.wcs.wcs.cdelt[0]
        beamsigma2 = self.header['BMIN'] / self.wcs.wcs.cdelt[0]
        return (np.pi * beamsigma1 * beamsigma2) / (4 * np.log(2))

