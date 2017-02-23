import logging
import math
import numpy as np


def matcher(catalog, sources, wcshelper):
    """
    We establish a source match if two conditions hold:
    1. the source ellipse overlaps the 1 sigma region of the catalog
       by more than 90%
    2. the source ellipse exceeds the 3 sigma region of the catalog
       by no more than 100%

    Returns:
        Either the index of the matching source or throws Exception
        if there is no suitable match.
    """
    # Set up grid to test match against 1 sigma overlap
    xmin, xmax, ymin, ymax = catalog.bounds()
    X, Y = np.mgrid[xmin:xmax+1, ymin:ymax+1]

    catalog_05sigma = catalog.draw(X, Y, magnify=0.5)
    count = np.sum(catalog_05sigma)

    # Best match
    best = None

    for i, source in enumerate(sources):
        # For speed, we first test for loose locality:
        # somewhere in the viscinity of 3 sigma of the catalog source
        max_distance = max(catalog.major, catalog.minor) * 8
        distance = math.sqrt((catalog.loc_x - source.loc_x)**2 + (catalog.loc_y - source.loc_y)**2)
        if distance > max_distance:
            continue

        # It's  close enough, so we go ahead and draw pixels and test for
        # sufficient overlap.
        source_box = source.draw(X, Y)
        missed = np.sum(catalog_05sigma[-source_box])
        if missed/count < 0.05:
            # Now we test for significant overreach beyond 4 sigma
            # First, we create a new grid that is the larger of the two:
            # either the source or the catalog at 4 sigma width
            catalog_bounds = catalog.bounds(magnify=4)
            source_bounds = source.bounds()
            xmin = min(catalog_bounds[0], source_bounds[0])
            xmax = max(catalog_bounds[1], source_bounds[1])
            ymin = min(catalog_bounds[2], source_bounds[2])
            ymax = max(catalog_bounds[3], source_bounds[3])
            X4, Y4 = np.mgrid[xmin:xmax+1, ymin:ymax+1]

            catalog_4sigma = catalog.draw(X4, Y4, magnify=4)
            count4 = np.sum(catalog_4sigma)
            source_box = source.draw(X4, Y4)
            excess = np.sum(source_box[-catalog_4sigma])
            # For small sources, we consider a bounding box of up to 150px x 150px a match.
            if np.sum(source_box) > 22500 and excess/count4 > 1.0:
                logging.info("Match excluded: catalog {} and detection {}, but 4 sigma excess is {} of 4 sigma pixels"
                             .format(catalog.id, source.id,  excess/count4))
                if best:
                    logging.info("Clobbering previous close match!")
                best = i
            else:
                logging.info("Source matched: catalog {} and detection {}, with 0.5 sigma miss of {}, 4 sigma excess of {}"
                             .format(catalog.id,
                                     source.id,
                                     missed/count,
                                     excess/count4))
                return i, 'match'
        elif missed/count < 0.99:
            logging.info("Close match: catalog {} and detection {}, missing {} of 0.5 sigma pixels"
                         .format(catalog.id, source.id, missed/count))
            if best:
                logging.info("Clobbering previous close match!")
            best = i

    else:
        if best:
            return best, 'close'
        else:
            raise FailedMatchException("No match")


class FailedMatchException(Exception):
    pass
