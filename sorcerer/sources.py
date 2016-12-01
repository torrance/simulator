import csv
import math
from collections import namedtuple
import numpy as np


Source = namedtuple('Source', ['id', 'loc_x', 'loc_y', 'peak', 'major',
                               'minor', 'pa', 'ra', 'dec', 'wmajor',
                               'wminor', 'wpa', 'total', 'totalerr'])


def synthetic_source_generator(count, xsize, ysize, wcshelper):
    for i in range(count):
        loc_x = np.random.randint(0, high=xsize)
        loc_y = np.random.randint(0, high=ysize)
        peak = 1.5 + np.abs(np.random.normal(0, 1))
        major = np.abs(np.random.normal(0, max(0.02 * xsize, 35))) + 20
        minor = major * np.random.uniform(low=0.5, high=1.5)
        pa = np.random.uniform(low=0.0, high=360.0)
        ra, dec, wmajor, wminor, wpa = wcshelper.pix2sky_ellipse(
            (loc_x, loc_y), major, minor, pa
        )
        total = (2 * np.pi * peak * major * minor) / wcshelper.beamarea_pix()

        yield Source._make([
            i,       # ID
            loc_x,   # loc_x
            loc_y,   # loc_y
            peak,    # peak
            major,   # major
            minor,   # minor
            pa,      # pa
            ra,      # ra
            dec,     # dec
            wmajor,  # wmajor
            wminor,  # wminor
            wpa,     # wpa
            total,   # total
            0,       # totalerr
        ])


def is_overlapping_source(candidate, sources):
    """
    Return True if candidate is within some threshold viscinity
    of another source.
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
    r = max(candidate.major, candidate.minor) * 1.5
    xmin = candidate.loc_x - r
    xmax = candidate.loc_x + r
    ymin = candidate.loc_y - r
    ymax = candidate.loc_y + r
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
                sources.append(Source._make([
                    int(row[0]),          # ID
                    ellipse[0],           # loc_x
                    ellipse[1],           # loc_y
                    float(row[10]),       # peak
                    ellipse[2],           # major
                    ellipse[3],           # minor
                    ellipse[4],           # angle
                    float(row[6]),        # ra
                    float(row[8]),        # dec
                    float(row[14])/3600,  # wmajor
                    float(row[16])/3600,  # wminor
                    float(row[18]),       # wpa
                    float(row[12]),       # total
                    float(row[13]),       # totalerr
                ]))
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
            sources.append(Source._make(source))
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
                sources.append(Source._make([
                    int(words[0]),        # ID
                    ellipse[0],           # loc_x
                    ellipse[1],           # loc_y
                    float(words[22]),     # peak
                    ellipse[2],           # major
                    ellipse[3],           # minor
                    float(words[12]),     # pa
                    float(words[7]),      # ra
                    float(words[8]),      # dec
                    float(words[11])/60,  # wmajor
                    float(words[10])/60,  # wminor
                    wpa,                  # wpa
                    float(words[19]),     # total
                    float(words[20]),     # totalerr
                ]))
            except:
                print("Failed to parse line: {}".format(line))
                raise
    return sources
