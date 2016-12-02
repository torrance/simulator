from collections import namedtuple


Ellipse = namedtuple('Ellipse', ['ra', 'dec', 'wmajor', 'wminor', 'pa'])


class Annotation:
    """
    A simple class to create KVIS annotion files.
    Usage:
       with Annotation(filename) as annotator:
           annotator.write_ellipse(ellipse1)
    """
    def __init__(self, filename, append=False):
        if append:
            self.file = open(filename, 'a')
        else:
            self.file = open(filename, 'w')
            self.file.write("PA STANDARD\n")
            self.file.write("COORD W\n")
            self.file.write("COLOR BLUE\n")
            self.file.write("\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def write_ellipse(self, ellipse, label=None, comment=None):
        line = "ELLIPSE {} {} {} {} {} # {}\n"
        self.file.write(line.format(ellipse.ra, ellipse.dec, ellipse.wmajor, ellipse.wminor, ellipse.pa, comment))
        if label:
            line = "TEXT {} {} {} # {}\n"
            self.file.write(line.format(ellipse.ra, ellipse.dec, label, comment))
        self.file.write("\n")
