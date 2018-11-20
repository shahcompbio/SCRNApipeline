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
import pandas

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
            sce_experiment = SingleCellExperiment.fromRData(sce)
            _genes = set(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
            genes = genes.intersection(_genes)
        return list(genes)

    @staticmethod
    def load(rdata, assay, symbol, filter_cells, filter_transcripts, common_genes):
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        assert symbol in sce_experiment.rowData.keys()
        genes = list(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
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
        assert matrix_o.shape[0] == len(barcodes), "not the same"
        _barcodes = []
        for barcode, row in zip(barcodes,matrix_o):
            if list(row).count(0.0) != len(row):
                matrix_f.append(row)
                _barcodes.append(barcode)
        barcodes = _barcodes
        matrix_t = numpy.transpose(matrix_f)
        return matrix_t, barcodes

    @staticmethod
    def run_em(rdata, filename, prefix, additional, assay="counts", symbol="Symbol", filter_cells=True, filter_transcripts=True):
        # sce_experiment = SingleCellExperiment.fromRData(rdata)
        # assert symbol in sce_experiment.rowData.keys()
        # genes = list(map(lambda x: x.upper(), list(sce_experiment.rowData["Symbol"])))
        # barcodes = list(sce_experiment.colData["Barcode"])
        # matrix = sce_experiment.assays[assay]
        # rho = eval(open(rho_matrix, "r").read())
        # rho = GeneMarkerMatrix(rho)
        # def subset(gene,row):
        #     return gene in rho.genes
        # assert len(genes) == matrix.shape[0]
        # rows = matrix[numpy.array([subset(gene,row) for gene,row in zip(genes,matrix)])]
        # genes = set(genes).intersection(set(rho.genes))
        # matrix = numpy.transpose(rows.toarray())
        # if filter_cells:
        #     _matrix = []
        #     for row in matrix:
        #         if list(map(int, row)).count(0) > 0:
        #             _matrix.append(list(row))
        #     matrix = numpy.array(_matrix)
        #     _matrix = []
        # matrix_t = numpy.transpose(matrix)
        # if filter_transcripts:
        #     _genes = []
        #     for gene, col in zip(genes,matrix_t):
        #         if list(map(int, row)).count(0) > 0:
        #             _matrix.append(list(col))
        #             _genes.append(gene)
        #     genes = _genes
        #     matrix = numpy.array(_matrix)
        # matrix_o = numpy.transpose(matrix)
        # matrix_f = []
        # assert matrix_o.shape[0] == len(barcodes), "not the same"
        # _barcodes = []
        # for barcode, row in zip(barcodes,matrix_o):
        #     if list(row).count(0.0) != len(row):
        #         matrix_f.append(row)
        #         _barcodes.append(barcode)
        # barcodes = _barcodes
        # matrix_t = numpy.transpose(matrix_f)
        rho = eval(open(rho_matrix, "r").read())
        rho = GeneMarkerMatrix(rho)
        all_rdata = [rdata] + additional
        genes = CellAssign.common_genes(all_rdata, symbol)
        genes = set(rho.genes).intersection(set(genes))
        matrix_t, barcodes = CellAssign.load(rdata, assay, symbol, filter_cells, filter_transcripts, genes)
        print("First", matrix_t.shape)
        if additional is not None:
            for sce in additional:
                _matrix_t, _barcodes = CellAssign.load(sce, assay, symbol, filter_cells, filter_transcripts, genes)
                print("Other",_matrix_t.shape)
                matrix_t = numpy.hstack((matrix_t,_matrix_t))
                barcodes += _barcodes
        matrix_t = numpy.array(matrix_t)
        print("SHAPE",matrix_t.shape)
        s = EdgeRInterface.calcNormFactors(matrix_t, method="TMM")
        s = pandas2ri.ri2py(s)
        rho_binary_matrix = numpy.array(rho.matrix(subset=genes,include_other=False))

        # df_dict = dict()
        # for gene, row in zip(genes,matrix_t):
        #     df_dict[gene] = row
        # df = pandas.DataFrame(df_dict)
        # cor = df.corr()
        #
        # cor.loc[:,:] =  numpy.tril(cor.values, k=-1)
        # cor = cor.stack()
        # groups = cor[cor > 0.5]
        # print(groups)

        matrix = numpy.transpose(matrix_t)
        assert matrix.shape[1] == rho_binary_matrix.shape[0], "Dimensions between rho and expression matrix do not match!"
        fit = CellAssignInterface.cellassign_em(matrix, rho_binary_matrix, s=s, data_type="RNAseq")
        #pickle.dump(fit, open(filename,"wb"))
        #fit = pickle.load(open(filename,"rb"))
        pyfit = dict(zip(fit.names, list(fit)))
        pyfit["Barcode"] = barcodes
        # mle_params = dict(zip(pyfit["mle_params"].names, list(pyfit["mle_params"])))
        conversion = dict(zip(sorted(list(set(pyfit["cell_type"]))),rho.celltypes()))
        cells = []
        for assignment in list(pyfit["cell_type"]):
            cells.append(conversion[assignment])
        pyfit["cell_type"] = cells
        pickle.dump(pyfit, open(filename,"wb"))
        #celltypes(cells,"cell_types.png",[cell.split("_")[0] for cell in rho.cells])
        robjects.r.assign("cell_assign_fit", fit)
        robjects.r("saveRDS(cell_assign_fit, file='{}_cellassign.rdata')".format(prefix))
