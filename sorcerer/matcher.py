import logging
import numpy as np


def matcher(catalog, sources, wcshelper):
    """
    We establish a source match if two conditions hold:
    1. the source ellipse overlaps the 1 sigma region of the catalog
       by more than 95%
    2. the source ellipse exceeds the 3 sigma region of the catalog
       by no more than 80%

    Returns:
        Either the index of the matching source or throws Exception
        if there is no suitable match.
    """
    # Set up grid to test match against 1 sigma overlap
    xmin, xmax, ymin, ymax = catalog.bounds()
    X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]

    catalog_1sigma = catalog.draw(X, Y)
    count = np.sum(catalog_1sigma)

    for i, source in enumerate(sources):
        source_box = source.draw(X, Y)
        missed = np.sum(catalog_1sigma[-source_box])
        if missed/count < 0.05:
            # Now we test for significant overreach beyond 3 sigma
            # First, we create a new grid that is the larger of the two:
            # either the source or the catalog at 3 sigma width
            catalog_bounds = catalog.bounds(magnify=3)
            source_bounds = source.bounds()
            xmin = min(catalog_bounds[0], source_bounds[0])
            xmax = max(catalog_bounds[1], source_bounds[1])
            ymin = min(catalog_bounds[2], source_bounds[2])
            ymax = max(catalog_bounds[3], source_bounds[3])
            X3, Y3 = np.mgrid[xmin:xmax+1, ymin:ymax+1]

            catalog_3sigma = catalog.draw(X3, Y3, magnify=3)
            count3 = np.sum(catalog_3sigma)
            source_box = source.draw(X3, Y3)
            excess = np.sum(source_box[-catalog_3sigma])
            if excess/count3 > 0.8:
                logging.info("Match excluded: catalog {} and detection {}, but 3 sigma excess is {} of 3 sigma pixels"
                             .format(catalog.id, source.id,  excess/count3))
            else:
                logging.info("Source matched: catalog {} and detection {}, with 1 sigma miss of {}, 3 sigma excess of {}"
                             .format(catalog.id,
                                     source.id,
                                     missed/count,
                                     excess/count3))
                return i
        elif missed/count < 0.5:
            logging.info("Close match: catalog {} and detection {}, missing {} of 1 sigma pixels"
                         .format(catalog.id, source.id, missed/count))

    else:
        raise FailedMatchException("No match")


class FailedMatchException(Exception):
    pass
