import os
import subprocess



class QCReport(object):

    def __init__(self, output):
        self.output = output

    def generate_script(self):
        script = os.path.join(self.output,"quality_control.R")
        output = open(script,"w")
        output.write("library(scater)\n")
        output.write("library(SingleCellExperiment)\n")
        output.write("library(EnsDb.Hsapiens.v86)\n")
        output.write("library(DropletUtils)\n")
        output.write("library(stringr)\n")
        output.write("library(scran)\n")
        output.write("library(Rtsne)")
        output.write("\n\n")

        output.write("args = commandArgs(trailingOnly=TRUE)\n\n")

        output.write("rdata <- readRDS(args[1])\n")
        output.write("sce <- as(rdata, 'SingleCellExperiment')\n")
        output.write("location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column='SEQNAME', keytype='GENEID')\n")
        output.write("rowData(sce)$CHR <- location\n")

        output.write("e.out <- emptyDrops(counts(sce))\n")
        output.write("sum(e.out$FDR <= 0.01, na.rm=TRUE)\n")
        output.write("sce <- sce[,which(e.out$FDR <= 0.01)]\n")

        output.write("""sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=='MT'), Ribo=which(str_detect(rowData(sce)$Symbol, '^RP(L|S)'))))\n""")

        output.write("png('umi_distribution.png')\n")
        output.write("hist(sce$log10_total_counts, breaks=20, col='grey80',xlab='Log-total UMI count')\n")
        output.write("dev.off()\n")

        output.write("png('features_distribution.png')\n")
        output.write("hist(sce$log10_total_features_by_counts, breaks=20, col='grey80',xlab='Log-total number of expressed features')\n")
        output.write("dev.off()\n")

        output.write("png('mito_distribution.png')\n")
        output.write("hist(sce$pct_counts_Mito, breaks=20, col='grey80',xlab='Proportion of reads in mitochondrial genes')\n")
        output.write("dev.off()\n")

        output.write("png('ribo_distribution.png')\n")
        output.write("hist(sce$pct_counts_Ribo, breaks=20, col='grey80',xlab='Proportion of reads in ribosomal genes')\n")
        output.write("dev.off()\n")


        output.write("high.mito <- isOutlier(sce$pct_counts_Mito, nmads=6, type='higher')\n")
        output.write("sce <- sce[,!high.mito]\n")

        output.write("keep_genes <- nexprs(sce, byrow=TRUE) >= 4\n")
        output.write("sce <- sce[keep_genes,]\n")
        output.write("sce <- normalize(sce)\n")
        output.write("sizeFactors(sce) <- librarySizeFactors(sce)\n")

        # output.write("png('average.png')\n")
        # output.write("ave <- calcAverage(sce)\n")
        # output.write("rowData(sce)$AveCount <- ave\n")
        # output.write("hist(ave, col='grey80')\n")
        # output.write('dev.off()\n')
        #
        # output.write("png('freqvmean.png')\n")
        # output.write("hist(ave, col='grey80')\n")
        # output.write("dev.off()\n")

        output.write("fit <- trendVar(sce, use.spikes=FALSE)\n")
        output.write("decomp <- decomposeVar(sce, fit)\n")
        output.write("top.hvgs <- order(decomp$bio, decreasing=TRUE)\n")

        """
        batch <- rep(c("1", "2"), each=100)
        design <- model.matrix(~batch)
        alt.fit2 <- trendVar(sce, design=design)
        alt.decomp2 <- decomposeVar(sce, alt.fit)
        """


        output.write("png('mean_var.png.png')\n")
        output.write("plot(decomp$mean, decomp$total, xlab='Mean log-expression', ylab='Variance')\n")
        output.write("dev.off()\n")


        # output.write("pca_data <- prcomp(t(log1p(assay(sce))))\n")
        # output.write("tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)\n")
        # output.write("reducedDims(sce) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)\n")

        # output.write("null.dist <- correlateNull(ncol(sce), iter=1e5)\n")
        # Only using the first 200 genes as a demonstration.
        # output.write("cor.pairs <- correlatePairs(sce, subset.row=top.hvgs[1:200], null.dist=null.dist)\n")
        # output.write("head(cor.pairs)\n")
        # output.write("write.csv(cor.pairs, file = 'correlated.csv', row.names=FALSE)\n")

        output.write("snn.gr <- buildSNNGraph(sce)\n")
        output.write("clusters <- igraph::cluster_walktrap(snn.gr)\n")
        output.write("sce$Cluster <- factor(clusters$membership)\n")

        output.write("print('Plotting by perplexity...')\n\n")

        output.write("png('tsne.png')\n")
        output.write("plotTSNE(sce, colour_by='Cluster', rerun=TRUE)\n")
        output.write("dev.off()\n")

        # output.write("png('highestExprs.png')\n")
        # output.write("plotHighestExprs(sce, exprs_values = 'counts')\n")
        # output.write("dev.off()\n")

        output.write("saveRDS(sce, file=args[2])\n")

        output.close()
        return script

    def generate_filtered_script(self):
        pass
