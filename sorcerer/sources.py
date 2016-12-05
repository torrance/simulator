import math
import numpy as np
from sorcerer.models import draw_ellipse, draw_rectangle
from sorcerer.output import Ellipse


class EllipticalSource:
    def __init__(self, **kwargs):
        self.id = kwargs['ID']
        self.loc_x = kwargs['loc_x']
        self.loc_y = kwargs['loc_y']
        self.peak = kwargs['peak']
        self.major = kwargs['major']
        self.minor = kwargs['minor']
        self.pa = kwargs['pa']
        self.ra = kwargs['ra']
        self.dec = kwargs['dec']
        self.wmajor = kwargs['wmajor']
        self.wminor = kwargs['wminor']
        self.wpa = kwargs['wpa']
        self.total = kwargs['total']
        self.totalerr = kwargs['totalerr']

    def __getitem__(self, key):
        return (
            self.id,
            self.loc_x,
            self.loc_y,
            self.peak,
            self.major,
            self.minor,
            self.pa,
            self.ra,
            self.dec,
            self.wmajor,
            self.wminor,
            self.wpa,
            self.total,
            self.totalerr,
        )[key]

    def columns(self):
        return [
            'ID',
            'loc_x',
            'loc_y',
            'peak',
            'major',
            'minor',
            'pa',
            'ra',
            'dec',
            'wmajor',
            'wminor',
            'wpa',
            'total',
            'totalerr',
        ]

    def annotate(self, ann, contours=None):
        ann.write_text(self.ra, self.dec, self.id)
        ann.write_ellipse(self, comment=self.id)

        if contours:
            for i in np.linspace(*contours):
                ann.write_ellipse(
                    Ellipse(self.ra, self.dec, self.wmajor*i, self.wminor*i, self.pa),
                    comment=self.id
                )

    def draw(self, X, Y, magnify=1):
        return draw_ellipse(X, Y, self.loc_x, self.loc_y, self.major*magnify, self.minor*magnify, self.pa)

    def show(self, fig, color=None):
        """Draw self onto APLpy FITSFigure."""
        fig.show_ellipses(
            self.ra,
            self.dec,
            self.wminor,
            self.wmajor,
            angle=-self.wpa,
            edgecolor=color,
            linewidth=0.6
        )

    def bounds(self, magnify=1):
        r = max(self.major, self.minor) * magnify
        xmin = int(self.loc_x-r)
        xmax = int(self.loc_x+r)
        ymin = int(self.loc_y-r)
        ymax = int(self.loc_y+r)
        return xmin, xmax, ymin, ymax

    def mean_width(self):
        return (self.major + self.minor)/2


class RectangularSource:
    def __init__(self, **kwargs):
        self.id = kwargs['ID']
        self.loc_x = kwargs['loc_x']
        self.loc_y = kwargs['loc_y']
        self.peak = kwargs['peak']
        self.xmin = kwargs['xmin']
        self.xmax = kwargs['xmax']
        self.ymin = kwargs['ymin']
        self.ymax = kwargs['ymax']
        self.ra = kwargs['ra']
        self.dec = kwargs['dec']
        self.wtopleft = kwargs['wtopleft']
        self.wtopright = kwargs['wtopright']
        self.wbottomright = kwargs['wbottomright']
        self.wbottomleft = kwargs['wbottomleft']
        self.total = kwargs['total']
        self.totalerr = kwargs['totalerr']

    def __getitem__(self, key):
        return (
            self.id,
            self.loc_x,
            self.loc_y,
            self.peak,
            self.xmin,
            self.xmax,
            self.ymin,
            self.ymax,
            self.ra,
            self.dec,
            self.wtopleft,
            self.wtopright,
            self.wbottomright,
            self.wbottomleft,
            self.total,
            self.totalerr,
        )[key]

    def columns(self):
        return [
            'ID',
            'loc_x',
            'loc_y',
            'peak',
            'xmin',
            'xmax',
            'ymin',
            'ymax',
            'ra',
            'dec',
            'wtopleft',
            'wtopright',
            'wbottomright',
            'wbottomleft',
            'total',
            'totalerr',
        ]

    def annotate(self, ann):
        ann.write_text(self.ra, self.dec, self.id)
        ann.write_clines([self.wtopleft, self.wtopright, self.wbottomright, self.wbottomleft, self.wtopleft])

    def draw(self, X, Y):
        return draw_rectangle(X, Y, self.xmin, self.xmax, self.ymin, self.ymax)

    def show(self, fig, color=None):
        """Draw self onto APLpy FITSFigure."""
        fig.show_polygons(
            [np.array([self.wtopleft, self.wtopright, self.wbottomright, self.wbottomleft])],
            edgecolor=color,
            linewidth=0.6
        )

    def bounds(self):
        return self.xmin, self.xmax, self.ymin, self.ymax

    def mean_width(self):
        return ((self.xmax - self.xmin)/2 + (self.ymax - self.ymin))/2


def synthetic_source_generator(count, xsize, ysize, wcshelper):
    for i in range(count):
        loc_x = np.random.randint(0, high=xsize)
        loc_y = np.random.randint(0, high=ysize)
        peak = np.random.uniform(low=1.5, high=5)
        major = np.random.uniform(low=20, high=200)
        minor = major * np.random.uniform(low=0.5, high=1.5)
        pa = np.random.uniform(low=0.0, high=360.0)
        ra, dec, wmajor, wminor, wpa = wcshelper.pix2sky_ellipse(
            (loc_x, loc_y), major, minor, pa
        )
        total = (2 * np.pi * peak * major * minor) / wcshelper.beamarea_pix()

        yield EllipticalSource(
            ID=i,
            loc_x=loc_x,
            loc_y=loc_y,
            peak=peak,
            major=major,
            minor=minor,
            pa=pa,
            ra=ra,
            dec=dec,
            wmajor=wmajor,
            wminor=wminor,
            wpa=wpa,
            total=total,
            totalerr=0,
        )


def is_overlapping_source(candidate, sources):
    """
    Return True if candidate is within some threshold viscinity
    of another source.

    Args:
        candidate: An EllipticalSource instance
        sources: A list of EllipticalSource instances
    """
    for source in sources:
        min_distance = (max(candidate.major, candidate.minor)
                        + max(source.major, source.minor)) * 3
        distance = math.sqrt((candidate.loc_x - source.loc_x)**2
                             + (candidate.loc_y - source.loc_y)**2)
        if distance < min_distance:
            return True
    return False


def is_edge(candidate, xsize, ysize):
    xmin, xmax, ymin, ymax = candidate.bounds(magnify=2)
    return xmin < 0 or xmax > xsize or ymin < 0 or ymax > ysize
