import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from interface.singlecellexperiment import SingleCellExperiment
from interface.genemarkermatrix import GeneMarkerMatrix

CellAssignInterface = importr("cellassign")
BiobaseInterface    = importr("Biobase")
SummarizedExperimentInterface = importr('SummarizedExperiment')
SingleCellExperimentInterface = importr('SingleCellExperiment')

class CellAssign(object):

    def __init__(self):
        pass

    def run_em(self, sce_experiment):
        rho = GeneMarkerMatrix()
        rho_binary_matrix = rho.load_example()
        sme = sce_experiment.asSummarizedExperiment()
        print(CellAssignInterface.cellassign_em(sme, rho_binary_matrix))
