library(clonealign)
library(SummarizedExperiment)
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="matrix file name", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$file)
read.matrix(opt$file, header = FALSE, sep = "\t", skip = 0)

# data(example_sce)
# N <- ncol(example_sce)
# G <- nrow(example_sce)
# C <- 3
#
# L <- rowData(example_sce)[, c("A", "B", "C")]
# print(L)
#
# clone_align_fit <- clonealign(example_sce, L, max_iter_em = 5)
#
# print(clone_align_fit)
