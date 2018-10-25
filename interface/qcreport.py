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
        output.write("library(stringr)")
        output.write("\n\n")

        output.write("args = commandArgs(trailingOnly=TRUE)\n\n")

        output.write("sce <- readRDS(args[1])\n")
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

        output.write("png('average.png')\n")
        output.write("ave <- calcAverage(sce)\n")
        output.write("rowData(sce)$AveCount <- ave\n")
        output.write("hist(ave, col='grey80')\n")
        output.write('dev.off()\n')

        output.write("png('freqvmean.png')\n")
        output.write("hist(ave, col='grey80')\n")
        output.write("dev.off()\n")

        output.write("png('highestExprs.png')\n")
        output.write("plotHighestExprs(sce, exprs_values = 'counts')\n")
        output.write("dev.off()\n")

        output.write("saveRDS(sce, file=args[2])\n")

        output.close()
        return script

    def generate_filtered_script(self):
        pass


    # def run(self):
    #     script = self._generate_script()
    #     print(script)
    #     subprocess.call(["Rscript", script])
