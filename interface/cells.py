
import collections

class Cell(object):

    def __init__(self):
        self.type = None
        self.barcode = None
        self.data = collections.defaultdict(float)


class Cells(object):

    def __init__(self):
        self.cells = list()

    def add_cell(self, cell):
        self.cells.append(cell)
