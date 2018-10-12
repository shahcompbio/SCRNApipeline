from software.dropletutils import DropletUtils
from interface.singlecellexperiment import SingleCellExperiment
from software.scater import Scater
from interface.scateranalysis import ScaterAnalysis
from utils import config

# import rpy2.robjects as robjects
# from rpy2.robjects.packages import importr
# from rpy2.robjects import r, pandas2ri
# from rpy2.rinterface import RRuntimeError
#
# pandas2ri.activate()
#
# ScaterInterface = importr("scater")


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
    def barcodeRanks(sce):
        utils = DropletUtils()
        ranks = utils.barcodeRanks(sce.assays["counts"])
        return dict(zip(ranks.names, map(list,list(ranks))))

    @staticmethod
    def emptyDrops(counts):
        utils = DropletUtils()
        cells = utils.emptyDrops(counts)
        cells = cells.slots["listData"]
        return dict(zip(cells.names, map(list,list(cells))))

    @staticmethod
    def calculateCPM(counts):
        scater = Scater()
        return scater.calculateCPM(counts)

    @staticmethod
    def calculateTPM(counts):
        scater = Scater()
        return scater.calculateCPM(counts)


    @staticmethod
    def calculateFPKM(counts):
        scater = Scater()
        return scater.calculateCPM(counts)


    @staticmethod
    def calculateQCMetrics(sce):
        scater = Scater()
        return scater.calculateQCMetrics(sce)

    @staticmethod
    def librarySizeFactors(counts):
        scater = Scater()
        return scater.librarySizeFactors(counts)

    @staticmethod
    def normalizeExprs(sce):
        scater = Scater()
        return scater.normalizeExprs(sce)

    @staticmethod
    def analysis(sce):
        output = "/Users/ceglian/Development/scater_output"
        sa = ScaterAnalysis(output)
        counts = sce.assays["counts"]
        exit(0)
        bcrank = TenX.barcodeRanks(sce)
        sa.barcodeRanks(bcrank)
        qcdsce = TenX.calculateQCMetrics(sce)
        print("RowData")
        for att in qcdsce.rowData.keys():
            print(att)
        print("ColData")
        for att in qcdsce.colData.keys():
            print(att)
