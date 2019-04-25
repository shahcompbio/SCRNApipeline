import subprocess
import numpy
import os
import sys
import pickle

"""
SINGULARITYENV_NPY_MKL_FORCE_INTEL=GNU SINGULARITYENV_PATH=$PATH:/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin/ SINGULARITYENV_R_HOME=/usr/local/lib/R/ SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/R/lib"""

""""
KLRF1 TRAC VIM
"""

from interface.genemarkermatrix import GeneMarkerMatrix

class CellAssign(object):

    @staticmethod
    def cmd(rdata, rho_csv, results, lsf=True):
        CellAssign.script(rdata, rho_csv, results)
        env = os.environ.copy()
        env["NPY_MKL_FORCE_INTEL"] = "GNU"
        env["PATH"] = os.environ["PATH"] + ":/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin/"
        env["R_HOME"] = "/usr/local/lib/R/"
        env["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":/usr/local/lib/R/lib"
        env["RETICULATE_PYTHON"] = "/home/ceglian/anaconda/bin/python3"
        submit = ["/home/ceglian/anaconda/bin/Rscript",".cache/run_cellassign.R"]
        #subprocess.call(submit, env=env)
        matched_results = os.path.join(os.path.split(results)[0],"cell_types.tsv")
        submit = ["/home/ceglian/anaconda/bin/Rscript",".cache/match.R"]
        subprocess.call(submit, env=env)

    @staticmethod
    def run(rdata, rho_yaml, results, rho_csv=".cache/rho.csv", lsf=True):
        if not os.path.exists(".cache"):
            os.makedirs(".cache")
        marker_list = GeneMarkerMatrix.read_yaml(rho_yaml)
        marker_list.write_matrix(rho_csv)
        CellAssign.cmd(rdata, rho_csv, results, lsf=lsf)
        print ("CellAssign finished.")
        matched_results = os.path.join(os.path.split(results)[0],"cell_types.tsv")
        pkl_fit = os.path.join(os.path.split(results)[0],"cell_types.pkl")
        lines = open(matched_results,"r").read().splitlines()
        header = lines.pop(0)
        barcodes = []
        celltypes = []
        pyfit = dict()
        for line in lines:
            line = [x.replace('"','') for x in line.split(",")]
            barcodes.append(line[1])
            celltypes.append(line[2])
        pyfit["Barcode"] = barcodes
        pyfit["cell_type"] = celltypes
        pickle.dump(pyfit, open(pkl_fit,"wb"))
        print ("Results written.")

    @staticmethod
    def script(rdata, rho_csv, results):
        filtered_sce = os.path.join(os.path.split(rdata)[0],"sce_cas.rdata")
        filtered_rho = os.path.join(os.path.split(rdata)[0],"rho_cas.rdata")
        matched_results = os.path.join(os.path.split(results)[0],"cell_types.tsv")
        configured = open(".cache/run_cellassign.R","w")
        configured.write(script.format(sce=rdata,rho=rho_csv,fname=results,fsce=filtered_sce,frho=filtered_rho))
        configured.close()
        match = open(".cache/match.R","w")
        match.write(match_barcodes.format(sce=filtered_sce,fit=results,fname=matched_results))
        match.close()


script = """
library(reticulate)
use_python("/home/ceglian/anaconda/bin/python3")
library(cellassign)
library(tensorflow)
library(scran)

rho <- read.csv("{rho}")
rownames(rho) <- rho$X
rho <- rho[,-1]

sce <- readRDS("{sce}")

# cells_to_keep <- sce$pct_counts_mito < 15
# table_cells_to_keep <- table(cells_to_keep)
# sce <- sce[,cells_to_keep]
# summ <- summary(sce$total_counts)
# thresh <- summ[[2]]
# keep <- sce$total_counts > thresh
# sce <- sce[,keep]

rownames(sce) <- rowData(sce)$Symbol
rho <- as.matrix(rho)
counts(sce) <- data.matrix(counts(sce))
sce <- sce[rowSums(counts(sce)) > 0,]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
sce <- sce[common_genes,]
sce <- sce[,colSums(counts(sce))>0]
sce <- sce[max(counts(sce)) > 1,]
sce <- sce[rowSums(counts(sce)) > 1,]
common_genes <- intersect(rowData(sce)$Symbol,rownames(rho))
rho <- rho[common_genes,]

rho <- data.matrix(rho)
s <- sizeFactors(sce)

library(tensorflow)

fit_cellassign <- cellassign(exprs_obj = sce, marker_gene_info = rho, s = s, shrinkage=TRUE)
saveRDS(fit_cellassign, file = '{fname}')
saveRDS(sce, file="{fsce}")
saveRDS(sce, file="{frho}")

"""


match_barcodes = """
library(SingleCellExperiment)
library(cellassign)
sce = readRDS("{sce}")
fit = readRDS("{fit}")

barcodes <- data.frame(barcode=colData(sce)$Barcode,celltype=fit$cell_type)
write.csv(barcodes, file="{fname}")

"""


if __name__ == '__main__':
    rho_yaml = "/work/shah/reference/transcriptomes/markers/hgsc_v1.yaml"
    rdata = "/work/shah/ceglian/Project_09443_D/ABDOM-CD45N_IGO_09443_D_2/runs/.cache/ABDOM-CD45N_IGO_09443_D_2/ABDOM-CD45N_IGO_09443_D_2.rdata"
    results = "cellassignfit.rds"
    CellAssign.run(rdata, rho_yaml, results)
