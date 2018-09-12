import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from singlecellexperiment import SingleCellExperiment


class DropletUtils(object):

    """

    """

    def __init__(self):
        self.namespace = importr("DropletUtils")
        self.tools = importr("scater")
        for attr, reference in self.namespace.__dict__.items():
            setattr(self, attr, reference)




# def read10xcounts(path_to_mtx):
#     tenx = importr("DropletUtils")
#     for f in tenx.__dict__.keys():
#         print(f)
#     res = tenx.read10xCounts(path_to_mtx)
#
#
# def annotaterows(singlecellexperiment):
#     utils = importr("DropletUtils")
#     print(utils.__dict__)

if __name__ == '__main__':
    tenx = DropletUtils()
    res = tenx.read10xCounts("~/data/raw_gene_bc_matrices/GRCh38/")
    sce = SingleCellExperiment.fromRS4(res)
    print(sce.rowData())
    print(sce.colData())
    print(sce.assayData().keys())
    print(sce.reducedDims())
    print(sce.sizeFactors())
    sce = SingleCellExperiment.fromRData("~/data/example_sce.RData")
    print(sce.assayData()["BatchCellMeans"].shape)
