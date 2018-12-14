import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.rinterface import RRuntimeError
from interface.singlecellexperiment import SingleCellExperiment
from interface.genemarkermatrix import GeneMarkerMatrix
from interface.tenxanalysis import TenxAnalysis
from utils.plotting import celltypes
import subprocess
import numpy
import os
import sys
from utils.config import *
import pickle
import pandas
import scanpy.api as sc

os.environ["RETICULATE_PYTHON"] = sys.executable

CellAssignInterface = importr("cellassign")
BiobaseInterface    = importr("Biobase")
SummarizedExperimentInterface = importr('SummarizedExperiment')
SingleCellExperimentInterface = importr('SingleCellExperiment')
EdgeRInterface = importr("edgeR")
ScaterInterface = importr("scater")

class CellAssign(object):

    @staticmethod
    def common_genes(sces, symbol):
        sce_experiment = SingleCellExperiment.fromRData(sces[0])
        genes = set(map(lambda x: x.upper(), list(sce_experiment.rowData[symbol])))
        for sce in sces[1:]:
            print(sce)
            sce_experiment = SingleCellExperiment.fromRData(sce)
            _genes = set(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
            genes = genes.intersection(_genes)
        return list(genes)

    @staticmethod
    def load(tenx, rdata, assay, symbol, filter_cells, filter_transcripts, common_genes):
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        assert symbol in sce_experiment.rowData.keys()
        genes = list(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
        convert = tenx.gene_map(genes)
        barcodes = list(sce_experiment.colData["Barcode"])
        matrix = sce_experiment.assays[assay]
        def subset(gene,row):
            return gene in common_genes
        assert len(genes) == matrix.shape[0]
        rows = matrix[numpy.array([subset(gene,row) for gene,row in zip(genes,matrix)])]
        genes = set(genes).intersection(set(common_genes))
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
        _barcodes = []
        for barcode, row in zip(barcodes,matrix_o):
            if list(row).count(0.0) != len(row):
                matrix_f.append(row)
                _barcodes.append(barcode)
        barcodes = _barcodes
        matrix_t = numpy.transpose(matrix_f)
        return matrix_t, barcodes, genes

    @staticmethod
    def write_input(matrix, rho, factors):
        robjects.r.assign("input_matrix",matrix)
        robjects.r("saveRDS(input_matrix,file='rdata/input_matrix.rdata')")
        robjects.r.assign("rho",rho)
        robjects.r("saveRDS(rho,file='rdata/rho.rdata')")
        robjects.r.assign("s",factors)
        robjects.r("saveRDS(s,file='rdata/factors.rdata')")

    @staticmethod
    def run_command_line(matrix, rho, factors):
        CellAssign.write_input(matrix, rho, factors)
        CellAssign.command()
        subprocess.call(["Rscript","run_cellassign.R"])

    @staticmethod
    def run_em(tenx, rdata, filename, prefix, rho_matrix, additional, assay="counts", symbol="Symbol", filter_cells=True, filter_transcripts=True, run_cmd=True, run_process=False):
        tenx = TenxAnalysis(tenx)
        sce = SingleCellExperiment.fromRData(rdata)
        print("Running EM.")
        if additional is not None:
            all_rdata = [rdata] + additional
        else:
            all_rdata = [rdata]
        sce = SingleCellExperiment.fromRData(rdata)
        genes = set(tenx.get_genes(sce))
        for nrdata in all_rdata[1:]:
            _genes = tenx.get_genes(SingleCellExperiment.fromRData(rdata))
            genes = genes.intersection(set(_genes))
        genes = list(genes)
        rho = GeneMarkerMatrix(rho_matrix)
        genes = set(rho.genes).intersection(set(genes))
        matrix_t, barcodes, origgenes = CellAssign.load(rdata, assay, symbol, filter_cells, filter_transcripts, genes)
        print(origgenes)
        if additional is not None:
            for sce in additional:
                _matrix_t, _barcodes, xgenes = CellAssign.load(sce, assay, symbol, filter_cells, filter_transcripts, genes)
                print(set(origgenes).difference(set(xgenes)))
                print(_matrix_t.shape, matrix_t.shape)
                matrix_t = numpy.hstack((matrix_t,_matrix_t))
                barcodes += _barcodes
        matrix_t = numpy.array(matrix_t)
        print("Calculating Size Factors.")
        s = EdgeRInterface.calcNormFactors(matrix_t, method="TMM")
        s = pandas2ri.ri2py(s)
        print("Finished size factors.")
        rho_binary_matrix = numpy.array(rho.matrix(subset=genes,include_other=False))
        print("Building Rho.")
        matrix = numpy.transpose(matrix_t)
        assert matrix.shape[1] == rho_binary_matrix.shape[0], "Dimensions between rho and expression matrix do not match!"
        if run_process:
            fit = CellAssignInterface.cellassign_em(exprs_obj = matrix, rho = rho_binary_matrix, s = s)
            robjects.r.assign("cell_assign_fit", fit)
            robjects.r("saveRDS(cell_assign_fit, file='rdata/cell_assign_fit.rdata')".format(prefix))
        else:
            print("Calling CMD Cell Assign.")
            if not os.path.exists("rdata/cell_assign_fit.rdata"):
                CellAssign.run_command_line(matrix, rho_binary_matrix, s)
            fit = r.readRDS("rdata/cell_assign_fit.rdata")
        pyfit = dict(zip(fit.names, list(fit)))
        pyfit["Barcode"] = barcodes
        conversion = dict(zip(sorted(list(set(pyfit["cell_type"]))),rho.celltypes()))
        cells = []
        for assignment in list(pyfit["cell_type"]):
            cells.append(conversion[assignment])
        pyfit["cell_type"] = cells
        pickle.dump(pyfit, open(filename,"wb"))

    @staticmethod
    def command():
        script = open("run_cellassign.R","w")
        script.write("""
        library(cellassign)
        gene_expression_data <- readRDS("./rdata/input_matrix.rdata")
        rho <- readRDS("./rdata/rho.rdata")
        s <- readRDS("./rdata/factors.rdata")
        cas <- cellassign_em(exprs_obj = gene_expression_data, rho = rho, s = s)
        saveRDS(cas,"./rdata/cell_assign_fit.rdata")
        """)
        script.close()
