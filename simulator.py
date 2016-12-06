#!/usr/bin/env python

import argparse
import csv
import logging
import math
import time
import aplpy
import matplotlib
import numpy as np
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from matplotlib import pyplot as plt
from sorcerer.models import draw_gaussian
from sorcerer.output import Annotation
from sorcerer.sources import synthetic_source_generator, is_overlapping_source, is_edge
from sorcerer.wcs_helpers import WCSHelper

# Use Ubuntu because it has greek glpyhs included (ie. mu)
matplotlib.rcParams['font.family'] = 'Ubuntu'
matplotlib.rcParams['font.size'] = '15.0'

parser = argparse.ArgumentParser(description='Generate a FITS radio image file.')
parser.add_argument('--xsize', '-x', dest='xsize', default=2000, type=int)
parser.add_argument('--ysize', '-y', dest='ysize', default=2000, type=int)
parser.add_argument('--count', '-c', dest='count', default=50, type=int)
parser.add_argument('--dist', dest='distribution', default='normal', choices=['normal', 'uniform'])
parser.add_argument('--image', dest='image', action='store_true')
parser.add_argument('--embed-original', dest='embed', action='store_true')
parser.add_argument('--log', '-l', dest='log', default='WARNING')
args = parser.parse_args()

XSIZE = args.xsize
YSIZE = args.ysize
COUNT = args.count
DISTRIBUTION = args.distribution
LOGLEVEL = args.log
IMAGE = args.image
EMBED = args.embed

filename = "simulated-" + time.strftime("%Y-%b-%d-%H%M")

# Configure logging
numeric_level = getattr(logging, LOGLEVEL.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: {}'.format(LOGLEVEL))
logging.basicConfig(level=numeric_level)

# Setup grid
data = np.zeros((YSIZE, XSIZE))

# Locate us in some junk location
w = wcs.WCS(naxis=2)
w.wcs.crpix = [int(XSIZE/2), int(YSIZE/2)]
w.wcs.crval = [0, 0]
w.wcs.cdelt = [0.0005, 0.0005]
w.wcs.ctype = ["RA---ZEA", "DEC--ZEA"]
w.wcs.equinox = 2000.0

# Add junk header information
logging.info("Creating FITs file...")
hdu0 = fits.PrimaryHDU(data)
hdu0.header['EXTEND'] = True
hdu0.header['BMAJ'] = 2 * (1/60)
hdu0.header['BMIN'] = 2 * (1/60)
hdu0.header['BPA'] = 0
hdu0.header['BUNIT'] = 'JY/BEAM'
hdu0.header += w.to_header()
wcshelper = WCSHelper(hdu0.header)

# Create sources, and remove overlapping sources if requested.
sources = []
for i, candidate in enumerate(synthetic_source_generator(COUNT, XSIZE, YSIZE, wcshelper)):
    if not is_overlapping_source(candidate, sources) and not is_edge(candidate, XSIZE, YSIZE):
        sources.append(candidate)

logging.info("Started with {}, left with {} after removing overlaps and edge sources."
             .format(COUNT, len(sources)))

# Create a scatterplot of sources for early analysis
# x = size; y = intensity
logging.info("Creating scatter plot...")
plt.figure('metrics', figsize=(8, 4.5))
plt.subplot(1, 1, 1)

xs = [source.mean_width() for source in sources]
ys = [source.peak for source in sources]
plt.scatter(xs, ys, c='lightblue')

plt.xlabel('Mean size of 1 sigma (pixels)')
plt.ylabel('Intensity (sigma)')
plt.savefig(filename + '-scatter.pdf')

# Draw the sources onto the data array
for i, source in enumerate(sources):
    logging.info("Drawing source {}/{}...".format(i+1, len(sources)))
    draw_gaussian(source, data)

# Hold onto a copy of the data prior to adding noise
if EMBED:
    hdu1 = fits.ImageHDU(np.copy(data))
    hdu1.header += w.to_header()

# Add noise
# We want a sigma of 1 for background, but this will be significantly decreased after
# convolution by a known factor, which we add here.
logging.info("Applying noise...")
convolution_factor = 2 * 4 * math.sqrt(math.pi)
noise = np.random.normal(0, 1 * convolution_factor, (YSIZE, XSIZE))
data += noise

# Convolve against a 2D Gaussian
logging.info("Applying Gaussian convolution...")
kernel = Gaussian2DKernel(4)
data = convolve(data, kernel)

logging.info("Saving fits file...")
hdu0.data = data
if EMBED:
    hdulist = fits.HDUList([hdu0, hdu1])
else:
    hdulist = fits.HDUList([hdu0])
hdulist.writeto(filename + '.fits')

# Write out annotation files
with Annotation(filename + '.ann') as ann:
    for source in sources:
        source.annotate(ann, contours=[1.5, 3, 4])

# Write out catalog of Gaussian images
with open(filename + '.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(sources[0].columns())
    for source in sources:
        writer.writerow(source)

# Also create pretty image filename
if IMAGE:
    logging.info("Creating image file...")
    fig = aplpy.FITSFigure(hdulist[0])
    logging.info("Applying contours...")
    fig.show_contour(hdulist[1], levels=[1.01, 1.5, 2, 3, 5])
    logging.info("Generating colorscale image")
    fig.show_colorscale(cmap='inferno')
    fig.add_grid()
    logging.info("Done. Saving...")
    fig.save(filename + '.png', dpi=320)
