#!/usr/bin/env python

import argparse
import csv
import logging
import math
import time
import aplpy
import numpy as np
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from sorcerer.models import draw_gaussian, Ellipse
from sorcerer.output import Annotation
from sorcerer.sources import synthetic_source_generator, is_overlapping_source
from sorcerer.wcs_helpers import WCSHelper


parser = argparse.ArgumentParser(description='Generate a FITS radio image file.')
parser.add_argument('--xsize', '-x', dest='xsize', default=2000, type=int)
parser.add_argument('--ysize', '-y', dest='ysize', default=2000, type=int)
parser.add_argument('--count', '-c', dest='count', default=50, type=int)
parser.add_argument('--dist', dest='distribution', default='normal', choices=['normal', 'uniform'])
parser.add_argument('--dedup', '-d', dest='dedup', action='store_true')
parser.add_argument('--log', '-l', dest='log', default='WARNING')
args = parser.parse_args()

XSIZE = args.xsize
YSIZE = args.ysize
COUNT = args.count
DEDUP = args.dedup
DISTRIBUTION = args.distribution
LOGLEVEL = args.log

# Configure logging
numeric_level = getattr(logging, LOGLEVEL.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: {}'.format(LOGLEVEL))
logging.basicConfig(level=numeric_level)

# Setup grid
data = np.zeros((XSIZE, YSIZE))

# Locate us in some junk location
w = wcs.WCS(naxis=2)
w.wcs.crpix = [int(XSIZE/2), int(YSIZE/2)]
w.wcs.crval = [0, 0]
w.wcs.cdelt = [0.005, 0.005]
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
    if not DEDUP or not is_overlapping_source(candidate, sources):
        sources.append(candidate)

if DEDUP:
    logging.info("Started with {}, left with {} after deduping."
                 .format(COUNT, len(sources)))

# Draw the sources onto the data array
for i, source in enumerate(sources):
    logging.info("Drawing source {}/{}...".format(i+1, len(sources)))
    draw_gaussian(source, data)

# Hold onto a copy of the data prior to adding noise, so that we can use this for contours
hdu1 = fits.ImageHDU(np.copy(data))
hdu1.header += w.to_header()

# Add noise
# We want a sigma of 1 for background, but this will be significantly decreased after
# convolution by a known factor, which we add here.
logging.info("Applying noise...")
convolution_factor = 2 * 4 * math.sqrt(math.pi)
noise = np.random.normal(0, 1 * convolution_factor, (XSIZE, YSIZE))
data += noise

# Convolve against a 2D Gaussian
logging.info("Applying Gaussian convolution...")
kernel = Gaussian2DKernel(4)
data = convolve(data, kernel)

logging.info("Saving fits file...")
hdu0.data = data
hdulist = fits.HDUList([hdu0, hdu1])
filename = "simulated-" + time.strftime("%Y-%b-%d-%H%M")
hdulist.writeto(filename + '.fits')

# Write out annotation files
with Annotation(filename + '.ann') as ann:
    for source in sources:
        # Also write ellipse for 2 and 3 sigma boundary
        ann.write_ellipse(source, label=source.id, comment=source.id)
        for i in np.linspace(1.5, 3, 4):
            world = wcshelper.pix2sky_ellipse(
                (source.loc_x, source.loc_y),
                i*source.major,
                i*source.minor,
                source.pa
            )
            ann.write_ellipse(
                Ellipse(world[0], world[1], world[2], world[3], source.pa),
                comment=source.id
            )

# Write out catalog of Gaussian images
with open(filename + '.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['ID', 'Loc-x', 'Loc-y', 'Peak', 'Major', 'Minor',
                     'PA', 'RA', 'DEC', 'WMajor', 'WMinor', 'WPA',
                     'Total', 'TotalErr'])
    for source in sources:
        writer.writerow(source)

# Also create pretty image filename
logging.info("Creating image file...")
fig = aplpy.FITSFigure(hdulist[0])
logging.info("Applying contours...")
fig.show_contour(hdulist[1], levels=[1.01, 1.5, 2, 3, 5])
logging.info("Generating colorscale image")
fig.show_colorscale(cmap='inferno')
fig.add_grid()
logging.info("Done. Saving...")
fig.save(filename + '.png', dpi=320)
