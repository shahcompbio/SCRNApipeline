from software.dropletutils import DropletUtils
from interface.singlecellexperiment import SingleCellExperiment
from software.scater import Scater
from interface.scateranalysis import ScaterAnalysis
from interface.qcreport import QCReport
from utils import config
import os


class TenX(object):

    def __init__(self):
        pass

    @staticmethod
    def read10xCountsRaw(tenx_analysis, output):
        utils = DropletUtils()
        counts = utils.read10xCounts(tenx_analysis.raw_matrices())
        sce = SingleCellExperiment.fromRS4(counts)
        sce.save(output)
        return sce

    @staticmethod
    def read10xCountsFiltered(tenx_analysis, output):
        utils = DropletUtils()
        print("TENX PATH", tenx_analysis.filtered_matrices())
        counts = utils.read10xCounts(tenx_analysis.filtered_matrices())
        sce = SingleCellExperiment.fromRS4(counts)
        sce.save(output)
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
    def calculateQCMetrics(sce, prefix, output):
        scater = Scater()
        qc_sce = scater.calculateQCMetrics(sce)
        qc_sce.save(os.path.join(output,"sce_qc_{}.rdata".format(prefix)))

    @staticmethod
    def plotHighestExprs(sce, prefix, output):
        scater = Scater()
        plot = scater.plotHighestExprs(sce)
        png = os.path.join(output, "highestExprs_{}.png".format(prefix))
        scater.save(png, plot)

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
        output = "tests/scater_output/"
        sa = ScaterAnalysis(output)
        bcrank = TenX.barcodeRanks(sce)
        sa.barcodeRanks(bcrank)
        qcdsce = TenX.calculateQCMetrics(sce)
