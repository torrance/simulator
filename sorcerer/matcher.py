import logging
import numpy as np
from sorcerer.models import ellipse


def matcher(catalog, sources, wcshelper):
    """
    We establish a source match if two conditions hold:
    1. the source ellipse overlaps the 2 sigma region of the catalog
       by more than 98%
    2. the source ellipse exceeds the 3 sigma region of the catalog
       by no more than 50%

    Returns:
        Either the index of the matching source or False
    """
    # Set up grid to test match against 2 sigma overlap
    r = max(catalog.major, catalog.minor)*2
    xmin = int(catalog.loc_x-r)
    xmax = int(catalog.loc_x+r)
    ymin = int(catalog.loc_y-r)
    ymax = int(catalog.loc_y+r)
    X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]

    catalog_2sigma = ellipse(X, Y, catalog.loc_x, catalog.loc_y,
                             2*catalog.major, 2*catalog.minor, catalog.pa)
    count = np.sum(catalog_2sigma)

    for i, source in enumerate(sources):
        source_ellipse = ellipse(X, Y, source.loc_x, source.loc_y,
                                 source.major, source.minor, source.pa)
        missed = np.sum(catalog_2sigma[-source_ellipse])
        if missed / count < 0.02:
            # Now we test for significant overreach beyond 3 sigma
            # First, we create a new grid that is the larger of the two:
            # either the source or the catalog at 3 sigma width
            r = max(source.major, source.minor, 3*catalog.major, 3*catalog.minor)
            xmin = int(source.loc_x-r)
            xmax = int(source.loc_x+r)
            ymin = int(source.loc_y-r)
            ymax = int(source.loc_y+r)
            X3, Y3 = np.mgrid[xmin:xmax+1, ymin:ymax+1]

            catalog_3sigma = ellipse(X3, Y3, catalog.loc_x, catalog.loc_y,
                                     3*catalog.major, 3*catalog.minor, catalog.pa)
            count3 = np.sum(catalog_3sigma)
            source_ellipse = ellipse(X3, Y3, source.loc_x, source.loc_y,
                                     source.major, source.minor, source.pa)
            excess = np.sum(source_ellipse[-catalog_3sigma])
            if excess/count3 > 0.5:
                logging.info("Match excluded: catalog {} and detection {}, \
                              but 3 sigma excess is {} of 3 sigma pixels"
                             .format(catalog.id, source.id,  excess/count3))
            else:
                logging.info("Source matched: catalog {} and detection {}, \
                              with 2 sigma miss of {}, 3 sigma excess of {}"
                             .format(catalog.id,
                                     source.id,
                                     missed/count,
                                     excess/count3))
                return i
        elif missed/count < 0.5:
            logging.info("Close match: catalog {} and detection {}, \
                          missing {} of 2 sigma pixels"
                         .format(catalog.id, source.id, missed/count))

    else:
        return False