#!/usr/bin/env python

import argparse
import csv
import logging
import os.path
from astropy import wcs
from astropy.io import fits
import aplpy
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from sorcerer.output import Annotation
from sorcerer.input import *
from sorcerer.matcher import matcher, FailedMatchException
from sorcerer.models import draw_ellipse
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
        xmin, xmax, ymin, ymax = source.bounds()
        xmin = max(0, xmin)
        xmax = min(xmax, maxX)
        ymin = max(0, ymin)
        ymax = min(ymax, maxY)
        X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]
        X = X.ravel()
        Y = Y.ravel()
        Z = source.draw(X, Y)
        data1[Y, X] = Z
    for source in sources:
        xmin, xmax, ymin, ymax = source.bounds()
        xmin = max(0, xmin)
        xmax = min(xmax, maxX)
        ymin = max(0, ymin)
        ymax = min(ymax, maxY)
        X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]
        X = X.ravel()
        Y = Y.ravel()
        Z = source.draw(X, Y)
        data2[Y, X] = Z

    hdu0 = fits.PrimaryHDU(data1)
    hdu0.header = header
    hdu1 = fits.ImageHDU(data2)
    hdu1.header += wcs.WCS(header).to_header()
    hdulist = fits.HDUList([hdu0, hdu1])
    try:
        os.remove('debug_shapes.fits')
    except:
        pass
    hdulist.writeto('debug_shapes.fits')


parser = argparse.ArgumentParser(description='Compare source detections against catalog.')
parser.add_argument('--catalog', '-c', dest='catalog', required=True)
parser.add_argument('--fitsfile', '-f', dest='fitsfile', required=True)
parser.add_argument('--debug', dest='debug', action='store_true')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--aegean', dest='aegean')
group.add_argument('--duchamp', dest='duchamp')
group.add_argument('--blobcat', dest='blobcat')
group.add_argument('--oddity', dest='oddity')
group.add_argument('--selavy', dest='selavy')
parser.add_argument('--image', dest='image', action='store_true')
parser.add_argument('--log', '-l', dest='log', default='WARNING')
args = parser.parse_args()

CATALOG = args.catalog
AEGEAN = args.aegean
DUCHAMP = args.duchamp
BLOBCAT = args.blobcat
ODDITY = args.oddity
SELAVY = args.selavy
FITSFILE = args.fitsfile
LOGLEVEL = args.log
IMAGE = args.image
DEBUG = args.debug

# Configure logging
numeric_level = getattr(logging, LOGLEVEL.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: {}'.format(LOGLEVEL))
logging.basicConfig(level=numeric_level)

wcshelper = WCSHelper(fits.getheader(FITSFILE))
if DEBUG or IMAGE:
    hdulist = fits.open(FITSFILE, memmap=True)

catalog = synthetic_sources_from_CSV(CATALOG)
if AEGEAN:
    filename = AEGEAN
    sources = aegean_sources_from_CSV(filename, wcshelper)
elif DUCHAMP:
    filename = DUCHAMP
    sources = duchamp_sources_from_txt(filename, wcshelper)
elif BLOBCAT:
    filename = BLOBCAT
    sources = blobcat_sources_from_txt(filename, wcshelper)
elif ODDITY:
    filename = ODDITY
    sources = oddity_sources_from_csv(filename, wcshelper)
elif SELAVY:
    filename = SELAVY
    sources = selavy_sources_from_txt(filename, wcshelper)

filename = os.path.splitext(os.path.basename(filename))[0]

# Debugging
if DEBUG:
    shape = hdulist[0].data.shape
    header = hdulist[0].header
    debug_shapes(catalog, sources, header, shape)

logging.info("Matching up sources...")
matchings = dict()  # Keys: catalogID; Values: sourceID
close_matches = dict()
for i, source in enumerate(catalog):
    try:
        match, match_type = matcher(source, sources, wcshelper)
        if match_type == 'match':
            matchings[i] = match
        elif match_type == 'close':
            close_matches[i] = match
    except FailedMatchException:
        pass

matches = [(catalog[k], sources[v]) for k, v in matchings.items()]  # (catalog, source)
close_matches = [(catalog[k], sources[v]) for k, v in close_matches.items()]
ghosts = [source for i, source in enumerate(sources) if i not in matchings.values()]
misses = [source for i, source in enumerate(catalog) if i not in matchings.keys()]
logging.info("Sources matched.")

print("Matched: {} Missed: {} Ghost: {}".format(len(matches), len(misses), len(ghosts)))
print("False positive rate: {}%".format(len(ghosts)/len(matches)*100))

# Output CSV of sources
with open(filename + '-matched.csv', 'w') as f:
    writer = csv.writer(f)
    # Creating CSV headers will fail if there are no matches at all.
    # We just blow up in that case.
    writer.writerow(matches[0][0].columns() + matches[0][1].columns())
    emptyCatalog = [None] * len(matches[0][0].columns())
    emptySource = [None] * len(matches[0][1].columns())
    for catalog, source in matches:
        writer.writerow([*catalog, *source, 'MATCH'])
    for source in misses:
        writer.writerow([*source, *emptySource, 'MISS'])
    for source in ghosts:
        writer.writerow([*emptyCatalog, *source, 'GHOST'])

# Create annotation file
with Annotation(filename + '-matched.ann') as ann:
    ann.set_color('GREEN')
    for _, source in matches:
        source.annotate(ann)
    ann.set_color('ORANGE')
    for source in misses:
        source.annotate(ann)
    ann.set_color('yellow')
    for source in ghosts:
        source.annotate(ann)

# Create a scatterplot of matches/misses
# x = size; y = intensity
logging.info("Creating scatter plot...")
plt.figure('metrics', figsize=(8, 4.5))
plt.subplot(1, 1, 1)

xs = [source.mean_width() for source in misses]
ys = [source.peak for source in misses]
plt.scatter(xs, ys, c='orange')

xs = [source[0].mean_width() for source in matches]
ys = [source[0].peak for source in matches]
plt.scatter(xs, ys, c='green')

# xs = [source.mean_width() for source in ghosts]
# ys = [source.peak for source in ghosts]
# plt.scatter(xs, ys, c='lightblue')

plt.xlabel('Mean size of 1 sigma (pixels)')
plt.ylabel('Intensity (sigma)')
plt.savefig(filename + '-scatter.pdf')

if IMAGE:
    # Plot close matches
    for source, close in close_matches:
        # First find the overall bounding box
        xmin1, xmax1, ymin1, ymax1 = source.bounds()
        xmin2, xmax2, ymin2, ymax2 = close.bounds()

        xmin = min(xmin1, xmin2)
        xmax = max(xmax1, xmax2)
        ymin = min(ymin1, ymin2)
        ymax = max(ymax1, ymax2)

        # Enlarge the bounding box just slightly
        xwidth = xmax - xmin
        if xwidth > 300:
            xmin -= (xwidth * 1.2 - xwidth) / 2
            xmax += (xwidth * 1.2 - xwidth) / 2
        else:
            center = xmin + xwidth // 2
            xmin = center - 150
            xmax = center + 150
        yheight = ymax - ymin
        if yheight > 300:
            ymin -= (yheight * 1.2 - yheight) / 2
            ymax += (yheight * 1.2 - yheight) / 2
        else:
            center = ymin + yheight // 2
            ymin = center - 150
            ymax = center + 150

        maxX, maxY = hdulist[0].data.shape
        maxX, maxY = maxX - 1, maxY - 1
        xmin, xmax, ymin, ymax = int(max(0, xmin)), int(min(xmax, maxX)), int(max(0, ymin)), int(min(ymax, maxY))

        window = hdulist[0].data[ymin:ymax, xmin:xmax]
        window = fits.PrimaryHDU(window)
        window.header = hdulist[0].header.copy()
        window.header['CRPIX1'] = window.header['CRPIX1'] - xmin
        window.header['CRPIX2'] = window.header['CRPIX2'] - ymin

        fig = aplpy.FITSFigure(window)
        fig.show_colorscale(cmap='viridis')

        source.show(fig, color='green')
        close.show(fig, color='orange')

        fig.save(filename + '-closematch-ID' + str(source.id) + '-' + str(close.id) + '.png')

    close_sources = [source[1] for source in close_matches]
    for source in ghosts:
        if source in close_sources:
            continue

        xmin, xmax, ymin, ymax = source.bounds()

        # Enlarge the bounding box just slightly
        xwidth = xmax - xmin
        if xwidth > 300:
            xmin -= (xwidth * 1.2 - xwidth) / 2
            xmax += (xwidth * 1.2 - xwidth) / 2
        else:
            center = xmin + xwidth // 2
            xmin = center - 150
            xmax = center + 150
        yheight = ymax - ymin
        if yheight > 300:
            ymin -= (yheight * 1.2 - yheight) / 2
            ymax += (yheight * 1.2 - yheight) / 2
        else:
            center = ymin + yheight // 2
            ymin = center - 150
            ymax = center + 150

        maxX, maxY = hdulist[0].data.shape
        maxX, maxY = maxX - 1, maxY - 1
        xmin, xmax, ymin, ymax = int(max(0, xmin)), int(min(xmax, maxX)), int(max(0, ymin)), int(min(ymax, maxY))

        window = hdulist[0].data[ymin:ymax, xmin:xmax]
        window = fits.PrimaryHDU(window)
        window.header = hdulist[0].header.copy()
        window.header['CRPIX1'] = window.header['CRPIX1'] - xmin
        window.header['CRPIX2'] = window.header['CRPIX2'] - ymin

        fig = aplpy.FITSFigure(window)
        fig.show_colorscale(cmap='viridis')

        source.show(fig, color='orange')
        fig.save(filename + '-ghost-ID' + str(source.id) + '.png')

    # # Let's plot the matches and misses
    # logging.info("Creating matched image file...")
    # fig = aplpy.FITSFigure(hdulist[0])
    # fig.show_colorscale(cmap='inferno')

    # # Plot the misses
    # for source in misses:
    #     source.show(fig, color='orange')

    # # Plot the matches
    # for _, source in matches:
    #     source.show(fig, color='green')

    # # Plot the ghosts
    # for source in ghosts:
    #     source.show(fig, color='lightblue')

    # fig.add_grid()
    # logging.info("Saving matched image...")
    # fig.save(filename + '-matched.png', dpi=320)
