import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from interface.singlecellexperiment import SingleCellExperiment

pandas2ri.activate()

ScaterInterface = importr('scater')


class Scater(object):


    def __init__(self):
        self.namespace = importr("scater")
        for attr, reference in self.namespace.__dict__.items():
            if hasattr(self,attr):
                setattr(self, "_{}".format(attr), reference)
            else:
                setattr(self, attr, reference)

    def calculateFPKM(self, assay):
        res = self._calculateFPKM(assay.toarray())
        return pandas2ri.ri2py(res)

    def calculateCPM(self, assay):
        res = self._calculateCPM(assay.toarray())
        return pandas2ri.ri2py(res)

    def calculateTPM(self, assay):
        res = self._calculateTPM(assay.toarray())
        return pandas2ri.ri2py(res)

    def calculateQCMetrics(self, sce):
        res = self._calculateQCMetrics(sce)
        return SingleCellExperiment.fromRS4(res)

    def librarySizeFactors(self, assay):
        res = self._librarySizeFactors(assay)
        return pandas2ri.ri2py(res)

    def normalizeExprs(self, sce):
        res = self._normalizeExprs(sce)
        return SingleCellExperiment.fromRS4(res)

    def save(self, filename, plot):
        robjects.r.ggsave(filename=filename, plot=plot)
