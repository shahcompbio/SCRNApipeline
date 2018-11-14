import os
import subprocess



class QCReport(object):

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
        script = os.path.join(self.output,"quality_control.R")
        output = open(script,"w")

        # output.write("library(scater)\n")
        # output.write("library(SingleCellExperiment)\n")
        # output.write("library(EnsDb.Hsapiens.v86)\n")
        # output.write("library(DropletUtils)\n")
        # output.write("library(stringr)\n")
        # output.write("library(scran)\n")
        # output.write("library(Rtsne)\n")
        # output.write("\n\n")
        self.imports(output)

        output.write("args = commandArgs(trailingOnly=TRUE)\n\n")

        # output.write("rdata <- readRDS(args[1])\n")
        # output.write("sce <- as(rdata, 'SingleCellExperiment')\n")
        self.read(output,"args[1]")

        # output.write("location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column='SEQNAME', keytype='GENEID')\n")
        # output.write("rowData(sce)$CHR <- location\n")
        self.annotate(output)


        # output.write("e.out <- emptyDrops(counts(sce))\n")
        # output.write("sum(e.out$FDR <= 0.01, na.rm=TRUE)\n")
        # output.write("sce <- sce[,which(e.out$FDR <= 0.01)]\n")
        self.filter_empty_drops(output)

        #output.write("""sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=='MT'), Ribo=which(str_detect(rowData(sce)$Symbol, '^RP(L|S)'))))\n""")
        self.qc_metrics(output)


        output.write("png('umi_distribution.png')\n")
        #output.write("hist(sce$log10_total_counts, breaks=20, col='darkgoldenrod1',xlab='Log-total UMI count')\n")
        self.umi(output)
        output.write("dev.off()\n")

        output.write("png('features_distribution.png')\n")
        output.write("hist(sce$log10_total_features_by_counts, breaks=20, col='grey80',xlab='Log-total number of expressed features')\n")
        output.write("dev.off()\n")

        output.write("png('mito_distribution.png')\n")
        #output.write("hist(sce$pct_counts_Mito, breaks=20, col='darkolivegreen4',xlab='Proportion of reads in mitochondrial genes')\n")
        self.mito(output)
        output.write("dev.off()\n")

        output.write("png('ribo_distribution.png')\n")
        #output.write("hist(sce$pct_counts_Ribo, breaks=20, col='firebrick4',xlab='Proportion of reads in ribosomal genes')\n")
        self.ribo(output)
        output.write("dev.off()\n")


        #output.write("high.mito <- isOutlier(sce$pct_counts_Mito, nmads=6, type='higher')\n")
        #output.write("sce <- sce[,!high.mito]\n")
        self.filter_high_mito(output)

        #output.write("keep_genes <- nexprs(sce, byrow=TRUE) >= 4\n")
        #output.write("sce <- sce[keep_genes,]\n")
        self.filter_low_count_genes(output)

        #output.write("sce <- normalize(sce)\n")
        self.normalize(output)

        #output.write("sizeFactors(sce) <- librarySizeFactors(sce)\n")
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
        #output.write("rowData(sce)$SymbolUnique <- make.names(rowData(sce)$Symbol, unique=T)\n")
        #output.write("plotHighestExprs(sce, exprs_values = 'counts', feature_names_to_plot= 'SymbolUnique')\n")
        self.highest_exprs(output)
        output.write("dev.off()\n")

        output.write("saveRDS(sce, file=args[2])\n")

        output.close()
        return script
