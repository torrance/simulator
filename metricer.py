#!/usr/bin/env python

import argparse
import logging
import os.path
from astropy import wcs
from astropy.io import fits
import aplpy
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from sorcerer.models import ellipse
from sorcerer.matcher import matcher
from sorcerer.sources import aegean_sources_from_CSV, synthetic_sources_from_CSV, duchamp_sources_from_txt
from sorcerer.wcs_helpers import WCSHelper

# Use Ubuntu because it has greek glpyhs included (ie. mu)
matplotlib.rcParams['font.family'] = 'Ubuntu'
matplotlib.rcParams['font.size'] = '15.0'


def debug_shapes(catalog, sources, header, shape):
    """
    This is an ugly function that spits out a new
    FITs file with two layers: the 2 sigma ellipses of the catalog
    and the ellipses of the sources. This allows us to verify visually
    that the ellipses are oriented correctly and as we would expect.
    """
    maxY, maxX = shape
    data1 = np.zeros((maxY, maxX))
    data2 = np.zeros((maxY, maxX))
    maxY, maxX = maxY-1, maxX-1  # Allow for zero-indexing

    for source in catalog:
        r = max(source.major, source.minor)*2
        xmin = max(0, int(source.loc_x-r))
        xmax = min(maxX, int(source.loc_x+r))
        ymin = max(0, int(source.loc_y-r))
        ymax = min(maxX, int(source.loc_y+r))
        X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]
        X = X.ravel()
        Y = Y.ravel()

        Z = ellipse(X, Y, source.loc_x, source.loc_y, 2*source.major, 2*source.minor, source.pa)
        Xdash = X[Z]
        Ydash = Y[Z]
        data1[Ydash, Xdash] = 1
    for source in sources:
        r = max(source.major, source.minor)
        xmin = max(0, int(source.loc_x-r))
        xmax = min(maxX, int(source.loc_x+r))
        ymin = max(0, int(source.loc_y-r))
        ymax = min(maxX, int(source.loc_y+r))
        X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]
        X = X.ravel()
        Y = Y.ravel()

        Z = ellipse(X, Y, source.loc_x, source.loc_y, source.major, source.minor, source.pa)
        Xdash = X[Z]
        Ydash = Y[Z]
        data2[Ydash, Xdash] = 1

    hdu0 = fits.PrimaryHDU(data1)
    hdu0.header = header
    hdu1 = fits.ImageHDU(data2)
    hdu1.header += wcs.WCS(header).to_header()
    hdulist1 = fits.HDUList([hdu0, hdu1])
    try:
        os.remove('debug_shapes.fits')
    except:
        pass
    hdulist1.writeto('debug_shapes.fits')


parser = argparse.ArgumentParser(description='Compare source detections against catalog.')
parser.add_argument('--catalog', '-c', dest='catalog', required=True)
parser.add_argument('--fitsfile', '-f', dest='fitsfile', required=True)
parser.add_argument('--debug', dest='debug', action='store_true')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--aegean', dest='aegean')
group.add_argument('--duchamp', dest='duchamp')
parser.add_argument('--log', '-l', dest='log', default='WARNING')
args = parser.parse_args()

CATALOG = args.catalog
AEGEAN = args.aegean
DUCHAMP = args.duchamp
FITSFILE = args.fitsfile
LOGLEVEL = args.log
DEBUG = args.debug

# Configure logging
numeric_level = getattr(logging, LOGLEVEL.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: {}'.format(LOGLEVEL))
logging.basicConfig(level=numeric_level)

hdulist = fits.open(FITSFILE)
wcshelper = WCSHelper(hdulist[0].header)

catalog = synthetic_sources_from_CSV(CATALOG)
if AEGEAN:
    filename = AEGEAN
    sources = aegean_sources_from_CSV(filename, wcshelper)
elif DUCHAMP:
    filename = DUCHAMP
    sources = duchamp_sources_from_txt(filename, wcshelper)
filename = os.path.splitext(filename)[0]

# Debugging
if DEBUG:
    shape = hdulist[0].data.shape
    header = hdulist[0].header
    debug_shapes(catalog, sources, header, shape)

logging.info("Matching up sources...")
matchings = dict()  # Keys: catalogID; Values: aegeanID
for i, source in enumerate(catalog):
    match = matcher(source, sources, wcshelper)
    if match:
        matchings[i] = match

matches = [(catalog[k], sources[v]) for k, v in matchings.items()]
ghosts = [source for i, source in enumerate(sources) if i not in matchings.values()]
misses = [source for i, source in enumerate(catalog) if i not in matchings.keys()]
logging.info("Sources matched.")

print("Matched: {} Missed: {} Ghost: {}".format(len(matches), len(misses), len(ghosts)))

# Let's plot the matches and misses
logging.info("Creating matched image file...")
fig = aplpy.FITSFigure(hdulist[0])
fig.show_colorscale(cmap='inferno')

# Plot the misses
if len(misses):
    xw = [source.ra for source in misses]
    yw = [source.dec for source in misses]
    fig.show_markers(xw, yw, c='orange')

# Plot the matches
if len(matches):
    xw = [source[1].ra for source in matches]
    yw = [source[1].dec for source in matches]
    aw = [source[1].wminor for source in matches]
    bw = [source[1].wmajor for source in matches]
    pw = [-source[1].wpa for source in matches]
    fig.show_ellipses(xw, yw, aw, bw, angle=pw, edgecolor='green', linewidth=1)

# Plot the ghosts
if len(ghosts):
    xw = [source.ra for source in ghosts]
    yw = [source.dec for source in ghosts]
    aw = [source.wminor for source in ghosts]
    bw = [source.wmajor for source in ghosts]
    pw = [-source.wpa for source in ghosts]
    fig.show_ellipses(xw, yw, aw, bw, angle=pw, edgecolor='lightblue', linewidth=1)

fig.add_grid()
logging.info("Saving matched image...")
fig.save(filename + '-matched.png', dpi=320)

# Also create a scatterplot of matches/misses
# x = size; y = intensity
logging.info("Creating scatter plot...")
plt.figure('metrics', figsize=(8, 4.5))
plt.subplot(1, 1, 1)
xs = [(source.major + source.minor)/2 for source in misses]
ys = [source.intensity for source in misses]
plt.scatter(xs, ys, c='orange')
xs = [(source[0].major + source[0].minor)/2 for source in matches]
ys = [source[0].intensity for source in matches]
plt.scatter(xs, ys, c='green')
plt.xlabel('Mean size of 1 sigma (pixels)')
plt.ylabel('Intensity (sigma)')
plt.savefig(filename + '-scatter.pdf')
