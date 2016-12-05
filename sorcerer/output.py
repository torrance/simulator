from collections import namedtuple


Ellipse = namedtuple('Ellipse', ['ra', 'dec', 'wmajor', 'wminor', 'pa'])


class Annotation:
    """
    A simple class to create KVIS annotion files.
    Usage:
       with Annotation(filename) as annotator:
           annotator.write_ellipse(ellipse1)
    """
    def __init__(self, filename, color='BLUE', pa='STANDARD', append=False):
        if append:
            self.file = open(filename, 'a')
        else:
            self.file = open(filename, 'w')
            self.file.write("PA {}\n".format(pa))
            self.file.write("COORD W\n")
            self.file.write("COLOR {}\n".format(color))
            self.file.write("\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    def write_ellipse(self, ellipse, comment=None):
        line = "ELLIPSE {} {} {} {} {} # {}\n"
        if comment:
            line += " # {}".format(comment)
        line += "\n"
        self.file.write(line.format(ellipse.ra, ellipse.dec, ellipse.wmajor, ellipse.wminor, ellipse.pa, comment))

    def write_clines(self, vertices, label=None, comment=None):
        line = "CLINES "
        for vertex in vertices:
            line += " {} {}".format(*vertex)
        if comment:
            line += " # {}".format(comment)
        line += "\n"
        self.file.write(line)

    def write_text(self, ra, dec, text, comment=None):
        line = "TEXT {} {} {}".format(ra, dec, text)
        if comment:
            line += " # {}".format(comment)
        line += "\n"
        self.file.write(line)

    def set_color(self, color):
        self.file.write("COLOR {}\n\n".format(color))
