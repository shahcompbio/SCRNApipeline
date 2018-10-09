import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.rinterface import RRuntimeError
from interface.singlecellexperiment import SingleCellExperiment
from interface.genemarkermatrix import GeneMarkerMatrix
import numpy
import os


os.environ["RETICULATE_PYTHON"] = "/usr/local/bin/python3"

CellAssignInterface = importr("cellassign")
BiobaseInterface    = importr("Biobase")
SummarizedExperimentInterface = importr('SummarizedExperiment')
SingleCellExperimentInterface = importr('SingleCellExperiment')
EdgeRInterface = importr("edgeR")

class CellAssign(object):

    def __init__(self):
        pass

    def run_em(self, sce_experiment, rho):
        rho_binary_matrix = rho.matrix()
        def subset(gene,row):
            return gene in rho.genes

        genes = sce_experiment.rowData[0]
        genes = ["Gene{}".format(gene) for gene in list(genes)]
        matrix = sce_experiment.assays["counts"]
        genes = list(genes)
        assert len(genes) == matrix.shape[0]
        rows = matrix[numpy.array([subset(gene,row) for gene,row in zip(genes,matrix)])]
        sce_experiment.assays["counts"] = rows
        sme = sce_experiment.asSummarizedExperiment()
        matrix = numpy.transpose(rows.toarray())
        s = EdgeRInterface.calcNormFactors(matrix,method="TMM")
        _matrix = []
        for row in matrix:
            if sum(row) != 0:
                row = list(row)
                print(row)
                _matrix.append(row)
        matrix = numpy.array(_matrix)
        _matrix = []
        for col in numpy.transpose(matrix):
            if sum(col) != 0:
                _matrix.append(list(col))
        matrix = numpy.array(_matrix)
        matrix = numpy.transpose(matrix)

        assert matrix.shape[1] == rho_binary_matrix.shape[0], "Dimensions between rho and expression matrix do not match!"
        fit = CellAssignInterface.cellassign_em(matrix, rho_binary_matrix)
