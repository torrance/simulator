import csv
import math
from collections import namedtuple
import numpy as np
from sorcerer.models import draw_ellipse
from sorcerer.output import Ellipse


Source = namedtuple('Source', ['id', 'loc_x', 'loc_y', 'peak', 'major',
                               'minor', 'pa', 'ra', 'dec', 'wmajor',
                               'wminor', 'wpa', 'total', 'totalerr'])


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

    def annotate(self, ann):
        ann.write_ellipse(self, label=self.id, comment=self.id)
        for i in np.linspace(1.5, 3, 4):
            ann.write_ellipse(
                Ellipse(self.ra, self.dec, self.wmajor*i, self.wminor*i, self.pa),
                comment=self.id
            )

    def draw(self, X, Y, magnify=1):
        return draw_ellipse(X, Y, self.loc_x, self.loc_y, self.major*magnify, self.minor*magnify, self.pa)

    def bounds(self, magnify=1):
        r = max(self.major, self.minor) * magnify
        xmin = int(self.loc_x-r)
        xmax = int(self.loc_x+r)
        ymin = int(self.loc_y-r)
        ymax = int(self.loc_y+r)
        return xmin, xmax, ymin, ymax


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


def aegean_sources_from_CSV(filename, wcshelper):
    sources = []
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # For now, we are just taking a single source per island,
            # and ignoring summits
            if row[1] == '0':
                # Note that PA is Sky
                ellipse = wcshelper.sky2pix_ellipse(
                    (float(row[6]), float(row[8])),
                    float(row[14])/3600,
                    float(row[16])/3600,
                    float(row[18])
                )
                sources.append(EllipticalSource(
                    ID=int(row[0]),
                    loc_x=ellipse[0],
                    loc_y=ellipse[1],
                    peak=float(row[10]),
                    major=ellipse[2],
                    minor=ellipse[3],
                    angle=ellipse[4],
                    ra=float(row[6]),
                    dec=float(row[8]),
                    wmajor=float(row[14])/3600,
                    wminor=float(row[16])/3600,
                    wpa=float(row[18]),
                    total=float(row[12]),
                    totalerr=float(row[13]),
                ))
    return sources


def synthetic_sources_from_CSV(filename):
    sources = []
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        for i, row in enumerate(reader):
            # Skip first row
            if i == 0:
                continue
            # We gotta convert to float
            source = []
            for i, val in enumerate(row):
                if i == 0:
                    source.append(int(val))
                else:
                    source.append(float(val))
            sources.append(EllipticalSource(
                ID=source[0],
                loc_x=source[1],
                loc_y=source[2],
                peak=source[3],
                major=source[4],
                minor=source[5],
                pa=source[6],
                ra=source[7],
                dec=source[8],
                wmajor=source[9],
                wminor=source[10],
                wpa=source[11],
                total=source[12],
                totalerr=source[13],
            ))
    return sources


def duchamp_sources_from_txt(filename, wcshelper):
    sources = []
    with open(filename) as f:
        for line in f:
            try:
                if line[0] == '#':
                    continue
                words = line.split()
                # Note that PA is Standard
                # and major and minor appear to be reversed.
                ellipse = wcshelper.sky2pix_ellipse(
                    (float(words[7]), float(words[8])),
                    float(words[11])/60,
                    float(words[10])/60,
                    float(words[12]),
                )
                _, _, _, _, wpa = wcshelper.pix2sky_ellipse(
                    (ellipse[0], ellipse[1]),
                    ellipse[2],
                    ellipse[3],
                    float(words[12])
                )
                sources.append(EllipticalSource(
                    ID=int(words[0]),
                    loc_x=ellipse[0],
                    loc_y=ellipse[1],
                    peak=float(words[22]),
                    major=ellipse[2],
                    minor=ellipse[3],
                    pa=float(words[12]),
                    ra=float(words[7]),
                    dec=float(words[8]),
                    wmajor=float(words[11])/60,
                    wminor=float(words[10])/60,
                    wpa=wpa,
                    total=float(words[19]),
                    totalerr=float(words[20]),
                ))
            except:
                print("Failed to parse line: {}".format(line))
                raise
    return sources
