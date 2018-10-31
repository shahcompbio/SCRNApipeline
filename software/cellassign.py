import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.rinterface import RRuntimeError
from interface.singlecellexperiment import SingleCellExperiment
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes
import numpy
import os
import sys
from utils.config import *
import pickle

os.environ["RETICULATE_PYTHON"] = sys.executable

CellAssignInterface = importr("cellassign")
BiobaseInterface    = importr("Biobase")
SummarizedExperimentInterface = importr('SummarizedExperiment')
SingleCellExperimentInterface = importr('SingleCellExperiment')
EdgeRInterface = importr("edgeR")
ScaterInterface = importr("scater")

class CellAssign(object):

    @staticmethod
    def run_em(rdata, filename, prefix, assay="counts", symbol="Symbol", filter_cells=True, filter_transcripts=True):
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        assert symbol in sce_experiment.rowData.keys()
        genes = list(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
        matrix = sce_experiment.assays[assay]
        rho = eval(open(rho_matrix, "r").read())
        rho = GeneMarkerMatrix(rho)
        def subset(gene,row):
            return gene in rho.genes
        assert len(genes) == matrix.shape[0]
        rows = matrix[numpy.array([subset(gene,row) for gene,row in zip(genes,matrix)])]
        genes = set(genes).intersection(set(rho.genes))
        matrix = numpy.transpose(rows.toarray())
        if filter_cells:
            _matrix = []
            for row in matrix:
                if list(map(int, row)).count(0) > 0:
                    _matrix.append(list(row))
            matrix = numpy.array(_matrix)
            _matrix = []
        matrix_t = numpy.transpose(matrix)
        if filter_transcripts:
            _genes = []
            for gene, col in zip(genes,matrix_t):
                if list(map(int, row)).count(0) > 0:
                    _matrix.append(list(col))
                    _genes.append(gene)
            genes = _genes
            matrix = numpy.array(_matrix)
        matrix_o = numpy.transpose(matrix)
        matrix_f = []
        for row in matrix_o:
            if list(row).count(0.0) != len(row):
                matrix_f.append(row)
        matrix_t = numpy.transpose(matrix_f)
        s = EdgeRInterface.calcNormFactors(matrix_t, method="TMM")
        s = pandas2ri.ri2py(s)
        matrix = numpy.transpose(matrix_t)
        rho_binary_matrix = numpy.array(rho.matrix(subset=genes))
        assert matrix.shape[1] == rho_binary_matrix.shape[0], "Dimensions between rho and expression matrix do not match!"
        fit = CellAssignInterface.cellassign_em(matrix, rho_binary_matrix, s=s, sce_assay=assay, data_type="RNAseq")
        pickle.dump(fit, open(filename,"wb"))
        pyfit = dict(zip(fit.names, list(fit)))
        mle_params = dict(zip(pyfit["mle_params"].names, list(pyfit["mle_params"])))
        conversion = dict(zip(sorted(list(set(pyfit["cell_type"]))),rho.celltypes()))
        cells = []
        for assignment in list(pyfit["cell_type"]):
            cells.append(conversion[assignment])
        celltypes(cells,"{}_cell_types.png".format(prefix),[cell.split("_")[0] for cell in rho.cells])
        robjects.r.assign("cell_assign_fit", fit)
        robjects.r("saveRDS(cell_assign_fit, file='{}_cellassign.rdata')".format(prefix))
