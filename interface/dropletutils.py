import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from interface.singlecellexperiment import SingleCellExperiment


class DropletUtils(object):

    def __init__(self):
        self.namespace = importr("DropletUtils")
        self.tools = importr("scater")
        for attr, reference in self.namespace.__dict__.items():
            setattr(self, attr, reference)


if __name__ == '__main__':
    tenx = DropletUtils()
    res = tenx.read10xCounts("~/data/raw_gene_bc_matrices/GRCh38/")
    sce = SingleCellExperiment.fromRS4(res)
    sce = SingleCellExperiment.fromRData("~/data/example_sce.RData")
    print(sce.assays["BatchCellMeans"])
