import subprocess
import numpy
import os
import sys
import pickle

from interface.genemarkermatrix import GeneMarkerMatrix

class CellAssign(object):

    @staticmethod
    def cmd(rdata, rho_csv, results):
        CellAssign.script(rdata, rho_csv, results)
        subprocess.call(["Rscript",".cache/run_cellassign.R"])

    @staticmethod
    def run(rdata, rho_yaml, results, rho_csv=".cache/rho.csv"):
        if not os.path.exists(".cache"):
            os.makedirs(".cache")
        marker_list = GeneMarkerMatrix.read_yaml(rho_yaml)
        marker_list.write_matrix(rho_csv)
        if not os.path.exists(results):
            CellAssign.cmd(rdata, rho_csv, results)
        print ("CellAssign finished.")
        # fit = r.readRDS(results)
        # pyfit = dict(zip(fit.names, list(fit)))
        # pyfit["Barcode"] = barcodes
        # conversion = dict(zip(sorted(list(set(pyfit["cell_type"]))),rho.celltypes()))
        # cells = []
        # for assignment in list(pyfit["cell_type"]):
        #     cells.append(conversion[assignment])
        # pyfit["cell_type"] = cells
        # pickle.dump(pyfit, open(filename,"wb"))

    @staticmethod
    def script(rdata, rho_csv, results):
        configured = open(".cache/run_cellassign.R","w")
        configured.write(script.format(sce=rdata,rho=rho_csv,fname=results))
        configured.close()

script = """
library(cellassign)
library(tensorflow)
library(cellassign.utils)
library(scran)

rho <- read.csv("{rho}")
rownames(rho) <- rho$X
rho <- rho[,-1]

sce <- readRDS("{sce}")

cells_to_keep <- sce$pct_counts_mito < 10
table_cells_to_keep <- table(cells_to_keep)
sce <- sce[,cells_to_keep]
rownames(sce) <- rowData(sce)$Symbol
rho <- as.matrix(rho)
counts(sce) <- data.matrix(counts(sce))
sce <- sce[rowSums(counts(sce)) > 0,]
sce <- sce[,colSums(counts(sce))>0]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
sce <- sce[common_genes,]
rho <- rho[common_genes,]
sce <- sce[rowSums(counts(sce)) > 0,]
sce <- sce[,colSums(counts(sce))>0]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
sce <- sce[common_genes,]
rho <- rho[common_genes,]

rho <- data.matrix(rho)
s <- sizeFactors(sce)

fit_cellassign <- cellassign(exprs_obj = sce, marker_gene_info = rho, s = s, B = 20,shrinkage = TRUE, verbose = TRUE,rel_tol_em = 1e-5, num_runs = 1, min_delta = 2)

saveRDS(fit_cellassign, file = '{fname}')"""


if __name__ == '__main__':
    rho_yaml = "/work/shah/reference/markers/hgsc_v1.yaml"
    rdata = "/work/shah/ceglian/Project_09443_D/ABDOM-CD45N_IGO_09443_D_2/runs/.cache/ABDOM-CD45N_IGO_09443_D_2/ABDOM-CD45N_IGO_09443_D_2.rdata"
    results = "cellassignfit.rds"
    CellAssign.run(rdata, rho_yaml, results)
