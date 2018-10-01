import pandas

class CloneAlignFit(object):

    def __init__(self, clones=None, genes=None, ml_params=None):
        self.genes = genes
        self.clones = pandas.read_csv(clones, sep="\t",skiprows=1,header=None)
