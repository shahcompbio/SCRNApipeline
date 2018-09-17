import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from interface.singlecellexperiment import SingleCellExperiment


class DropletUtils(object):

    def __init__(self):
        self.namespace = importr("DropletUtils")
        for attr, reference in self.namespace.__dict__.items():
            setattr(self, attr, reference)

    def ranks(self, sparse_matrix):
        dcg_matrix = SingleCellExperiment.CSRtoDCG(sparse_matrix)
        results = self.barcodeRanks(dcg_matrix)
        return results

    def emptydrops(self, sparse_matrix):
        dcg_matrix = SingleCellExperiment.CSRtoDCG(sparse_matrix)
        results = self.emptyDrops(dcg_matrix)
        return results

if __name__ == '__main__':
    pass
