library(clonealign)
library(SummarizedExperiment)

matrix <- as.matrix(read.table("./software/scripts/tmp_assay.tsv", sep = "\t", skip = 0))
colnames(matrix) <- NULL

N <- ncol(matrix)
G <- nrow(matrix)

C <- 3
L <- rowData(example_sce)[, c("A", "B", "C")]
clone_align_fit <- clonealign(example_sce, L, max_iter_em = 5)
write.table(clone_align_fit$clone, file="./software/scripts/tmp_clone_alignments.tsv", sep="\t")
#write.table(clone_align_fit, file="./software/scripts/tmp_ml_params.tsv", sep="\t")
