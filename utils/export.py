import os
import subprocess
from markdown2 import Markdown
import shutil

from interface.singlecellexperiment import SingleCellExperiment
from interface.fastqdirectory import FastQDirectory
from utils.config import Configuration

config = Configuration()

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

def exportFinalize(results):
    results.finalize()

def exportMD(results):
    results.finalize()
    markdown = os.path.join(results.report_dir,"{}_analysis.md".format(config.prefix))
    output = open(markdown,"w")
    output.write("# {} Analysis\n".format(config.prefix.replace('_',' ').upper()))
    output.write("***\n")
    output.write("\n\n## Pipeline Output\n")
    output.write(" - [CellRanger Summary]({})\n".format(results.summary))
    output.write(" - SCE Object (TenX Raw): unfiltered_sce.rdata\n\n".format(results.filtered_sce))
    output.write(" - SCE Object (TenX Filtered): filtered_sce.rdata\n\n".format(results.filtered_sce))
    output.write(" - SCE Object (QC'd): {}\n\n".format(results.filtered_sce))
    output.write(" - CellAssign R Object: {}\n\n".format(results.raw))
    output.write(" - QC Workflow {} - Standard\n\n".format(config.qc_type))
    output.write(" - R QC Workflow: {}\n\n".format(results.script))

    output.write("***\n")
    output.write("## Quality Control\n")
    code = open(results.script,"r").read()

    output.write("## Analysis\n")
    for plot in results.plots:
        output.write("\n### {}\n".format(plot["header"]))
        output.write("![{}]({})\n".format(plot["desc"],plot["path"]))

    output.close()
    convertHTML(markdown, results.report_dir)


def convertHTML(markdown, report):
    output = open(os.path.join(report, "{}_results.html".format(config.prefix)),"w")
    markdowner = Markdown()
    html = markdowner.convert(open(markdown,"r").read())
    output.write(html)
    output.close()

def exportRMD(fastq, analysis, scater_workflow, prefix, sce, output_path):
    with_code=False
    rmd = os.path.join(output_path,"{}_analysis.rmd".format(prefix))
    output = open(rmd,"w")
    output.write("# {} Analysis\n".format(prefix.replace('_',' ').capitalize()))
    output.write("***\n")
    output.write("\n\n## Pipeline Output\n")
    output.write(" - [CellRanger Summary]({})\n".format(analysis.summary()))
    output.write(" - FastQ location: {}\n\n".format(os.fastq.path))

    output.write("***\n")
    output.write("## Quality Control\n")
    output.write("### FastQC\n")
    for sample, summary in fastq.qc_reports():
        output.write(" - [{}]({})\n".format(sample, summary))

    umi = os.path.join(output_path,"umi_distribution.png")
    mito = os.path.join(output_path,"mito_distribution.png")
    ribo = os.path.join(output_path, "ribo_distribution.png")
    freq = os.path.join(output_path, "highestExprs.png")

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

    # if with_code: codeblock(output, scater_workflow.highest_exprs)
    # output.write("\n### Highly Expressed Transcripts\n")
    # output.write("![]({})\n".format(freq))

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
        output.write("library(Seurat)\n")
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

    def filter_cells(self, output):
        output.write('sce <- FilterCells(object = sce, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))\n')

    def qc_metrics(self, output):
        output.write("""sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=='MT'), Ribo=which(str_detect(rowData(sce)$Symbol, '^RP(L|S)'))))\n""")

    def umi(self, output):
        output.write("hist(sce$log10_total_counts, breaks=20, col='darkgoldenrod1',xlab='Log-total UMI count')\n")

    def mito(self, output):
        output.write("hist(sce$pct_counts_Mito, breaks=20, col='darkolivegreen4',xlab='Proportion of reads in mitochondrial genes')\n")

    def ribo(self, output):
        output.write("hist(sce$pct_counts_Ribo, breaks=20, col='firebrick4',xlab='Proportion of reads in ribosomal genes')\n")

    def mito_percentage(self, output):
        output.write("mito.genes <- grep(pattern = '^MT-', x = rownames(x = sce@data), value = TRUE)\n")
        output.write("percent.mito <- Matrix::colSums(sce@raw.data[mito.genes, ])/Matrix::colSums(sce@raw.data)\n")
        output.write("sce <- AddMetaData(object = sce, metadata = percent.mito, col.name = 'percent.mito')\n")

    def get_symbols(self, output):
        output.write("symbols <- rowData(sce)$Symbol\n")


    def set_symbols(self, output):
        output.write("rowData(sce)$Symbol <- symbols\n")

    # def set_symbols(self, output):
    #     output.write("sce <- AddMetaData(object = sce, metadata = symbols, row.name = 'Symbol')\n")

    def violin_gene_mito_umi(self, output):
        output.write('VlnPlot(object = sce, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)\n')

    def normalize_seurat(self, output):
        output.write('sce <- NormalizeData(object = sce, normalization.method = "LogNormalize", scale.factor = 10000)\n')

    def gene_plot(self, output):
        output.write('par(mfrow = c(1, 2))\n')
        output.write('GenePlot(object = sce, gene1 = "nUMI", gene2 = "percent.mito")\n')
        output.write('GenePlot(object = sce, gene1 = "nUMI", gene2 = "nGene")\n')

    def filter_high_mito(self, output, stds=6):
        output.write("high.mito <- isOutlier(sce$pct_counts_Mito, nmads={}, type='higher')\n".format(stds))
        output.write("sce <- sce[,!high.mito]\n")

    def filter_low_count_genes(self, output, n_genes=4):
        output.write("keep_genes <- nexprs(sce, byrow=TRUE) >= {}\n".format(n_genes))
        output.write("sce <- sce[keep_genes,]\n")

    def find_highly_variable(self, output):
        output.write('sce <- FindVariableGenes(object = sce, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)\n')
        output.write('highly_variable_genes <- sce@var.genes\n')

    def select_highly_variable(self, output):
        output.write("sce <- sce[highly_variable_genes,]\n")

    def regress_out(self, output):
        output.write('sce <- ScaleData(object = sce, vars.to.regress = c("nUMI", "percent.mito"))\n')

    def sce_to_seurat(self, output):
        output.write('sce <- Convert(from = sce, to = "seurat")\n')

    def seurat_to_sce(self, output):
        output.write('sce <- Convert(from = sce, to = "sce")\n')

    def add_dim_names(self, output):
        output.write('rownames(sce) <- make.names(rowData(sce)$Symbol, unique=TRUE)\n')
        output.write('colnames(sce) <- colData(sce)$Barcode\n')

    def normalize(self, output):
        output.write("sce <- normalize(sce)\n")

    def calc_size_factors(self, output):
        output.write("sizeFactors(sce) <- librarySizeFactors(sce)\n")

    def highest_exprs(self, output, make_unique = False):
        if make_unique:
            output.write("rowData(sce)$SymbolUnique <- make.names(rowData(sce)$Symbol, unique=T)\n")
            output.write("plotHighestExprs(sce, exprs_values = 'counts', feature_names_to_plot= 'SymbolUnique')\n")
        else:
            output.write("plotHighestExprs(sce, exprs_values = 'counts', feature_names_to_plot= 'Symbol')\n")

    def mean_variance_trend(self, output):
        pass
        # output.write("new.trend <- makeTechTrend(x=sce)\n")
        # output.write("fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))\n")
        # output.write("plot(fit$mean, fit$var, pch=16)\n")
        # output.write("curve(fit$trend(x), col='dodgerblue', add=TRUE)\n")
        # output.write("curve(new.trend(x), col='red', add=TRUE)\n")

    def plot_qc(self,output, log=False):
        if not log:
            output.write("plotQC(sce, type = 'highest-expression', exprs_values = 'counts')\n")
        else:
            output.write("plotQC(sce, type = 'highest-expression', exprs_values = 'logcounts')\n")

    def generate_script(self):
        try:
            shutil.rmtree("./rdata/h5")
        except Exception:
            pass
        script = os.path.join(self.output,"scater_workflow.R")
        output = open(script,"w")

        self.imports(output)
        output.write("args = commandArgs(trailingOnly=TRUE)\n\n")
        self.read(output,"args[1]")
        self.annotate(output)
        self.qc_metrics(output)

        self.filter_low_count_genes(output,n_genes=config.low_counts_genes_threshold)
        self.add_dim_names(output)
        self.normalize(output)

        output.write("png('./figures/highestExprs.png')\n")
        self.highest_exprs(output,make_unique=True)
        output.write("dev.off()\n")

        self.sce_to_seurat(output)
        self.mito_percentage(output)

        # output.write("png('./figures/gene_mito_umi.png')\n")
        # self.gene_plot(output)
        # output.write("dev.off()\n")
        #
        # output.write("png('./figures/violin.png')\n")
        # self.violin_gene_mito_umi(output)
        # output.write("dev.off()\n")

        # self.filter_cells(output)
        #
        # self.normalize_seurat(output)

        # output.write("png('./figures/highly_variable_genes.png')\n")
        # self.find_highly_variable(output)
        # output.write("dev.off()\n")

        self.regress_out(output)

        self.seurat_to_sce(output)
        # self.select_highly_variable(output)
        self.normalize(output)
        self.filter_high_mito(output,stds=config.stds)


        output.write("png('./figures/umi_distribution.png')\n")
        self.umi(output)
        output.write("dev.off()\n")
        output.write("png('./figures/features_distribution.png')\n")
        output.write("hist(sce$log10_total_features_by_counts, breaks=20, col='grey80',xlab='Log-total number of expressed features')\n")
        output.write("dev.off()\n")
        output.write("png('./figures/mito_distribution.png')\n")
        self.mito(output)
        output.write("dev.off()\n")
        output.write("png('./figures/ribo_distribution.png')\n")
        self.ribo(output)
        output.write("dev.off()\n")

        output.write("png('./figures/mean_variance_trend.png')\n")
        self.mean_variance_trend(output)
        output.write("dev.off()\n")


        output.write("print('running pca...')\n")
        output.write("sce <- runPCA(sce)\n")
        output.write("png('./figures/pca.png')\n")
        output.write("plotPCA(sce)\n")
        output.write("dev.off()\n")

        output.write("print('running umap...')\n")
        output.write("sce <- runUMAP(sce)\n")
        output.write("png('./figures/umap.png')\n")
        output.write("plotUMAP(sce)\n")
        output.write("dev.off()\n")

        output.write("print('running tsne...')\n")
        output.write("sce <- runTSNE(sce)\n")
        output.write("png('./figures/tsne.png')\n")
        output.write("plotTSNE(sce)\n")
        output.write("dev.off()\n")

        output.write("sce <- runPCA(sce,ncomponents=200)\n")
        output.write("rowData(sce)$Symbol <- rowData(sce)$gene\n")

        output.write("saveRDS(sce, file=args[2])\n")

        output.close()

        return script



if __name__ == '__main__':
    pass
