import csv
import numpy as np
from sorcerer.sources import EllipticalSource, RectangularSource


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
                    pa=ellipse[4],
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
                line = line.strip()
                if line[0] == '#':
                    continue
                words = line.split()
                # Note that PA is Standard
                # and major and minor appear to be reversed.
                ellipse = wcshelper.sky2pix_ellipse(
                    (float(words[7]), float(words[8])),
                    float(words[11])/3600,
                    float(words[10])/3600,
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
                    wmajor=float(words[11])/3600,
                    wminor=float(words[10])/3600,
                    wpa=wpa,
                    total=float(words[20]) / wcshelper.beamarea_pix(),
                    totalerr=float(words[20]),
                ))
            except:
                print("Failed to parse line: {}".format(line))
                raise
    return sources


def blobcat_sources_from_txt(filename, wcshelper):
    sources = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line[0] == '#':
                continue
            words = line.split()

            xmin = int(words[18])
            xmax = int(words[19])
            ymin = int(words[20])
            ymax = int(words[21])
            wtopleft = wcshelper.pix2sky((xmin, ymax))
            wtopright = wcshelper.pix2sky((xmax, ymax))
            wbottomright = wcshelper.pix2sky((xmax, ymin))
            wbottomleft = wcshelper.pix2sky((xmin, ymin))

            sources.append(RectangularSource(
                ID=int(words[0]),
                loc_x=float(words[2]),
                loc_y=float(words[3]),
                peak=float(words[32]),
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                ra=float(words[4]),
                dec=float(words[5]),
                wtopleft=wtopleft,
                wtopright=wtopright,
                wbottomright=wbottomright,
                wbottomleft=wbottomleft,
                total=float(words[37]),
                totalerr=float(words[38]),
            ))

    return sources


def oddity_sources_from_csv(filename, wcshelper):
    sources = []
    with open(filename) as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        next(reader)  # Skip units
        for row in reader:
            x1, y1, x2, y2 = float(row[1]), float(row[2]), float(row[3]), float(row[4])
            center = (np.mean([x1, x2]), np.mean([y1, y2]))
            wcenter = wcshelper.pix2sky(center)
            wtopleft = wcshelper.pix2sky((x1, y2))
            wtopright = wcshelper.pix2sky((x2, y2))
            wbottomright = wcshelper.pix2sky((x2, y1))
            wbottomleft = wcshelper.pix2sky((x1, y1))

            sources.append(RectangularSource(
                ID=int(row[0]),
                loc_x=center[0],
                loc_y=center[1],
                peak=float(row[8]),
                xmin=x1,
                xmax=x2,
                ymin=y1,
                ymax=y2,
                ra=wcenter[0],
                dec=wcenter[1],
                wtopleft=wtopleft,
                wtopright=wtopright,
                wbottomright=wbottomright,
                wbottomleft=wbottomleft,
                total=float(row[6]),
                totalerr=0,
            ))

    return sources


def selavy_sources_from_txt(filename, wcshelper):
    sources = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line[0] == '#':
                continue

            words = line.split()

            ellipse = wcshelper.sky2pix_ellipse(
                (float(words[5]), float(words[6])),
                float(words[14])/3600/2,
                float(words[15])/3600/2,
                float(words[16]),
            )

            _, _, _, _, wpa = wcshelper.pix2sky_ellipse(
                (ellipse[0], ellipse[1]),
                ellipse[2],
                ellipse[3],
                float(words[16])
            )

            sources.append(EllipticalSource(
                ID=words[0],
                loc_x=ellipse[0],
                loc_y=ellipse[1],
                peak=float(words[10]),
                major=ellipse[3],
                minor=ellipse[2],
                pa=float(words[16]),
                ra=float(words[5]),
                dec=float(words[6]),
                wmajor=float(words[15])/3600/2,
                wminor=float(words[14])/3600/2,
                wpa=wpa,
                total=float(words[12])/1000,
                totalerr=0,
            ))

    return sources
