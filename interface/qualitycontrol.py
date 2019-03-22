from interface.tenxanalysis import TenxAnalysis

from azure.common.client_factory import get_client_from_cli_profile
from azure.mgmt.compute import ComputeManagementClient
from azure.mgmt.network import NetworkManagementClient
from azure.storage.blob import BlockBlobService, PublicAccess

import subprocess
import os
import uuid

class QualityControl(object):

    def __init__(self, tenx, sampleid):
        self.tenx = tenx
        self.sampleid = sampleid
        self.sce = os.path.join(tenx.path,"{0}.rdata".format(sampleid))
        self.cache = ".cache"
        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        self.script = os.path.join(self.cache,"qc.R")
        if not os.path.exists(self.script):
            output = open(self.script,"w")
            output.write(script)
            output.write(filter)
            output.write(save)
            output.close()
        self.storage_account = "scrnadata"
        self.container = "rdataraw{}".format(self.tenx.detected_version)
        self.block_blob_service = BlockBlobService(account_name='scrnadata', sas_token='?sv=2018-03-28&ss=bfqt&srt=sco&sp=rwdlacup&se=2021-03-19T02:52:48Z&st=2019-02-22T19:52:48Z&spr=https&sig=4oAGvIyqi9LPY89t21iWWp4XbbIgpmaldgcwWUOuf14%3D')

    def run(self, mito=10):
        mat = self.tenx.filtered_matrices()
        subprocess.call(["Rscript", self.script, mat, self.sce, str(mito)])

    def build(self):
        self.build = os.path.join(self.cache,"build.R")
        if not os.path.exists(self.build):
            output = open(self.script,"w")
            output.write(script)
            output.write(save)
            output.close()
        mat = self.tenx.filtered_matrices()
        subprocess.call(["Rscript", self.script, mat, self.sce])

    def sce(self):
        return SingleCellExperiment.fromRData(self.sce)

    def upload(self):
        self.block_blob_service.create_blob_from_path(self.container,"{0}.rdata".format(self.sampleid), self.sce)


script = """
library(scater)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(DropletUtils)
library(stringr)
library(scran)

args = commandArgs(trailingOnly=TRUE)

sce <- read10xCounts(args[1])

# Need to get rid of hg19_ in front
rownames(sce) <- rownames(rowData(sce)) <- rowData(sce)$ID <- gsub("hg19_", "", rownames(sce))
rowData(sce)$Symbol <- gsub("hg19_", "", rowData(sce)$Symbol)


rowData(sce)$ensembl_gene_id <- rownames(sce)

sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id",
attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene",
"start_position", "end_position", "chromosome_name"),
dataset = "hsapiens_gene_ensembl")


# Calculate size factors
sce_result <- tryCatch({scran::computeSumFactors(sce)},error=function(e) {NULL})
poolsize <- 100;
while (is.null(sce_result) && poolsize >= 0) {
    sce_result <- tryCatch({scran::computeSumFactors(sce,sizes=poolsize)},error=function(e) {NULL})
    if (is.null(sce_result)) {
      poolsize <- poolsize - 10;
    }
}

# Compute log normal expression values
sce <- normalize(sce)


# Get Mitochondrial genes for QC:
mt_genes <- which(rowData(sce)$chromosome_name == "MT")
ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                      ribo = rownames(sce)[ribo_genes])

# Calculate QC metrics
sce <- calculateQCMetrics(sce, feature_controls = feature_ctrls)

# Make sure colnames are unique
colnames(sce) <- paste0(metadata(sce)$id, "_", sce$Barcode)

sce <- runPCA(sce, ncomponents = 3)
set.seed(123L)
sce <- runTSNE(sce)

"""


filter = """
cells_to_keep <- sce$total_features_by_counts > 1000 & sce$pct_counts_mito < args[3]
table_cells_to_keep <- table(cells_to_keep)
sce <- sce[,cells_to_keep]
"""

save = """
saveRDS(sce, file=args[2])
"""
