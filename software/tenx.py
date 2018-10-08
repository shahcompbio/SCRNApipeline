from software.dropletutils import DropletUtils
from interface.singlecellexperiment import SingleCellExperiment
from utils import config

class TenX(object):

    def __init__(self):
        pass

    @staticmethod
    def read10xCounts(tenx_analysis):
        utils = DropletUtils()
        counts = utils.read10xCounts(tenx_analysis.raw_matrices(config.build))
        sce = SingleCellExperiment.fromRS4(counts)
        return sce

    @staticmethod
    def barcodeRanks(tenx_analysis):
        utils = DropletUtils()

    @staticmethod
    def emptyDrops(tenx_analysis):
        utils = DropletUtils()
