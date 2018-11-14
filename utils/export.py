import os
import subprocess

from interface.singlecellexperiment import SingleCellExperiment
from interface.fastqdirectory import FastQDirectory

def cat(fastq, output, prefix):
    result = os.path.join(output, "{}.fastq.gz".format(prefix))
    cmd = ["cat"]
    for fq in fastq.get_fastqs:
        cmd.append(fq)
    cmd.append(">")
    cmd.append(result)
    print(" ".join(cmd))
    subprocess.call(cmd)


def codeblock(handle, callback):
    handle.write("\n```\n")
    callback(handle)
    handle.write("```\n\n")

def imports(handle):
    handle.write("\n\n")
    handle.write("```\n")
    handle.write("library(scater)\n")
    handle.write("library(scran)\n")
    handle.write("library(SingleCellExperiment)\n")
    handle.write("library(Rtsne)\n")
    handle.write("```\n")


def exportRMD(fastq, analysis, scater_workflow, prefix, sce, output_path):
    with_code=False
    rmd = os.path.join(output_path,"{}_analysis.rmd".format(prefix))
    output = open(rmd,"w")
    output.write("# {} Analysis\n".format(prefix.replace('_',' ').capitalize()))
    output.write("***\n")
    output.write("\n\n## Pipeline Output\n")
    output.write(" - [CellRanger Summary]({})\n".format(analysis.summary()))
    output.write(" - FastQ location: {}\n\n".format(fastq.path))

    output.write("***\n")
    output.write("## Quality Control\n")
    output.write("### FastQC\n")
    for sample, summary in fastq.qc_reports():
        output.write(" - [{}]({})\n".format(sample, summary))

    umi = os.path.join(output_path,"umi_distribution.png")
    mito = os.path.join(output_path,"mito_distribution.png")
    ribo = os.path.join(output_path, "ribo_distribution.png")
    freq = os.path.join(output_path, "highestExprs.png")
    # tsne = os.path.join(output, "tsne.png")
    # umap = os.path.join(output, "umap.png")
    tsne_by_cluster = os.path.join(output_path, "tsne_by_cluster.png")
    tsne_by_cell_type = os.path.join(output_path, "tsne_by_cell_type.png")
    cell_type_by_cluster = os.path.join(output_path, "cluster_by_cell_type.png")
    celltypes = os.path.join(output_path, "cell_types.png")

    if with_code:
        codeblock(output, scater_workflow.imports)

        codeblock(output, lambda x : scater_workflow.read(x, sce))

        codeblock(output, scater_workflow.umi)
    output.write("\n### UMI\n")
    output.write("![]({})\n".format(umi))

    if with_code: codeblock(output, scater_workflow.ribo)
    output.write("\n### Ribosomal\n")
    output.write("![]({})\n".format(ribo))

    if with_code: codeblock(output, scater_workflow.mito)
    output.write("\n### Mitochondrial\n")
    output.write("![]({})\n".format(mito))

    if with_code: codeblock(output, scater_workflow.highest_exprs)
    output.write("\n### Highly Expressed Transcripts\n")
    output.write("![]({})\n".format(freq))

    output.write("\n***\n")
    output.write("\n### TSNE (Clusters)\n")
    output.write("![]({})\n".format(tsne_by_cluster))
    output.write("\n### TSNE (Cell Type)\n")
    output.write("![]({})\n".format(tsne_by_cell_type))
    output.write("\n### Cell Types by Cluster\n")
    output.write("![]({})\n".format(cell_type_by_cluster))
    output.write("\n### Cell Types\n")
    output.write("\n![]({})\n".format(celltypes))

    output.close()


class ScaterCode(object):

    def __init__(self, output):
        self.output = output

    def imports(self, output):
        output.write("library(scater)\n")
        output.write("library(SingleCellExperiment)\n")
        output.write("library(EnsDb.Hsapiens.v86)\n")
        output.write("library(DropletUtils)\n")
        output.write("library(stringr)\n")
        output.write("library(scran)\n")
        output.write("library(Rtsne)\n")
        output.write("\n\n")

    def read(self, output, rdata):
        output.write("rdata <- readRDS({})\n".format(rdata))
        output.write("sce <- as(rdata, 'SingleCellExperiment')\n")

    def annotate(self, output):
        output.write("location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column='SEQNAME', keytype='GENEID')\n")
        output.write("rowData(sce)$CHR <- location\n")

    def filter_empty_drops(self, output):
        output.write("e.out <- emptyDrops(counts(sce))\n")
        output.write("sum(e.out$FDR <= 0.01, na.rm=TRUE)\n")
        output.write("sce <- sce[,which(e.out$FDR <= 0.01)]\n")

    def qc_metrics(self, output):
        output.write("""sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=='MT'), Ribo=which(str_detect(rowData(sce)$Symbol, '^RP(L|S)'))))\n""")

    def umi(self, output):
        output.write("hist(sce$log10_total_counts, breaks=20, col='darkgoldenrod1',xlab='Log-total UMI count')\n")

    def mito(self, output):
        output.write("hist(sce$pct_counts_Mito, breaks=20, col='darkolivegreen4',xlab='Proportion of reads in mitochondrial genes')\n")

    def ribo(self, output):
        output.write("hist(sce$pct_counts_Ribo, breaks=20, col='firebrick4',xlab='Proportion of reads in ribosomal genes')\n")

    def filter_high_mito(self, output):
        output.write("high.mito <- isOutlier(sce$pct_counts_Mito, nmads=6, type='higher')\n")
        output.write("sce <- sce[,!high.mito]\n")

    def filter_low_count_genes(self, output):
        output.write("keep_genes <- nexprs(sce, byrow=TRUE) >= 4\n")
        output.write("sce <- sce[keep_genes,]\n")

    def normalize(self, output):
        output.write("sce <- normalize(sce)\n")

    def calc_size_factors(self, output):
        output.write("sizeFactors(sce) <- librarySizeFactors(sce)\n")

    def highest_exprs(self, output):
        output.write("rowData(sce)$SymbolUnique <- make.names(rowData(sce)$Symbol, unique=T)\n")
        output.write("plotHighestExprs(sce, exprs_values = 'counts', feature_names_to_plot= 'SymbolUnique')\n")

    def generate_script(self):
        script = os.path.join(self.output,"scater_workflow.R")
        output = open(script,"w")

        self.imports(output)
        output.write("args = commandArgs(trailingOnly=TRUE)\n\n")
        self.read(output,"args[1]")
        self.annotate(output)
        self.filter_empty_drops(output)
        self.qc_metrics(output)
        output.write("png('umi_distribution.png')\n")
        self.umi(output)
        output.write("dev.off()\n")
        output.write("png('features_distribution.png')\n")
        output.write("hist(sce$log10_total_features_by_counts, breaks=20, col='grey80',xlab='Log-total number of expressed features')\n")
        output.write("dev.off()\n")
        output.write("png('mito_distribution.png')\n")
        self.mito(output)
        output.write("dev.off()\n")
        output.write("png('ribo_distribution.png')\n")
        self.ribo(output)
        output.write("dev.off()\n")
        self.filter_high_mito(output)
        self.filter_low_count_genes(output)
        self.normalize(output)
        self.calc_size_factors(output)

        #
        # output.write("fit <- trendVar(sce, use.spikes=FALSE)\n")
        # output.write("decomp <- decomposeVar(sce, fit)\n")
        # output.write("top.hvgs <- order(decomp$bio, decreasing=TRUE)\n")
        #
        # """
        # batch <- rep(c("1", "2"), each=100)
        # design <- model.matrix(~batch)
        # alt.fit2 <- trendVar(sce, design=design)
        # alt.decomp2 <- decomposeVar(sce, alt.fit)
        # """
        #
        #
        # output.write("png('mean_var.png.png')\n")
        # output.write("plot(decomp$mean, decomp$total, xlab='Mean log-expression', ylab='Variance')\n")
        # output.write("dev.off()\n")


        # output.write("pca_data <- prcomp(t(log1p(assay(sce))))\n")
        # output.write("tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)\n")
        # output.write("reducedDims(sce) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)\n")

        # output.write("null.dist <- correlateNull(ncol(sce), iter=1e5)\n")
        # output.write("cor.pairs <- correlatePairs(sce, subset.row=top.hvgs[1:200], null.dist=null.dist)\n")
        # output.write("head(cor.pairs)\n")
        # output.write("write.csv(cor.pairs, file = 'correlated.csv', row.names=FALSE)\n")

        output.write("snn.gr <- buildSNNGraph(sce)\n")
        output.write("clusters <- igraph::cluster_walktrap(snn.gr)\n")
        output.write("sce$Cluster <- factor(clusters$membership)\n")

        output.write("sce <- runPCA(sce)\n")
        output.write("png('pca.png')\n")
        output.write("plotPCA(sce)\n")
        output.write("dev.off()\n")

        output.write("sce <- runTSNE(sce)\n")
        output.write("png('tsne.png')\n")
        output.write("plotTSNE(sce)\n")
        output.write("dev.off()\n")
        output.write("png('tsne_by_cluster.png')\n")
        output.write("plotTSNE(sce, colour_by='Cluster')\n")
        output.write("dev.off()\n")

        output.write("png('highestExprs.png')\n")
        output.write("rowData(sce)$SymbolUnique <- make.names(rowData(sce)$Symbol, unique=T)\n")
        output.write("plotHighestExprs(sce, exprs_values = 'counts', feature_names_to_plot= 'SymbolUnique')\n")
        self.highest_exprs(output)
        output.write("dev.off()\n")

        output.write("saveRDS(sce, file=args[2])\n")

        output.close()

        return script



if __name__ == '__main__':
    pass
