import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import numpy
import os
import sys
import pandas
import collections
from pybedtools import BedTool, Interval
import subprocess

from utils.config import Configuration

from interface.singlecellexperiment import SingleCellExperiment

# os.environ["RETICULATE_PYTHON"] = sys.executable
pandas2ri.activate()

config = Configuration()

class CloneAlign(object):

    @staticmethod
    def run(tenx, rdata, copy_number_data, clone_assignments, assay="counts", run_cmd=True):
        sce = SingleCellExperiment.fromRData(rdata)
        _genes = tenx.get_genes(sce)
        convert = tenx.gene_map(sce)
        genes = []
        for gene in _genes:
            if gene in convert:
                genes.append(convert[gene])
            else:
                genes.append(gene)
        adata = tenx.create_scanpy_adata(sce)
        matrix = adata.X.T
        assert assay in sce.assayNames, "Assay not present in SCE."
        if run_cmd:
            print("Calling CMD Clone Align.")
            if not os.path.exists("rdata/clone_align.rdata"):
                CloneAlign.run_command_line(adata, copy_number_data, clone_assignments, genes)
            if not os.path.exists("rdata/clone_align.rdata"):
                raise ValueError("Rscript 'run_clonealign.r' Failed.")
            cal = r.readRDS("rdata/cell_assign_fit.rdata")
        else:
            CloneAlignInterface = importr("clonealign")
            cal = CloneAlignInterface.clonealign(matrix, cnv_data)
            robjects.r.assign("clone_align_fit", clone_align_fit)
            robjects.r("saveRDS(clone_align_fit, file='{}')".format(filename))

    @staticmethod
    def assemble_copy_number_data(copy_number_data, clone_assignments, gene_subset):
        subset = lambda feature: Interval(feature.chrom,feature.start,feature.stop,strand=feature.strand,name=feature.name)
        rows = open(clone_assignments,"r").read().splitlines()
        header = rows.pop(0)
        index = collections.defaultdict(list)
        for row in rows:
            row = row.split(",")
            index[row[1]].append(row[0])
        copy_number_data = pandas.read_csv(copy_number_data, header=0)
        cnv_by_clone = pandas.DataFrame()
        for clone, cells in index.items():
            cnv_by_clone[clone] = copy_number_data[cells].sum(axis=1)
            copy_number_data.drop(cells,axis=1, inplace=True)
        copy_number_data.drop(["width"],axis=1,inplace=True)
        copy_number_data["name"] = ["Segment_{}".format(i) for i in copy_number_data.index]
        cnv_by_clone["Segment"]  = ["Segment_{}".format(i) for i in copy_number_data.index]
        cnv_by_clone.astype(str, copy=False)
        copy_number_data.to_csv("./tables/cnv.bed",sep="\t",index=False,header=False)
        cnv = BedTool("./tables/cnv.bed")
        genes = BedTool(config.genes_gtf)
        genes = genes.each(subset)
        annotated = genes.intersect(cnv, wb=True).moveto("./tables/annotated.bed")
        annotation_df = pandas.read_csv("./tables/annotated.bed",sep="\t",header=None,dtype=str)
        annotation_df.columns = ["GeneChrom","GeneStart","GeneEnd","GeneSymbol","GeneScore","GeneStrand","SegmentChrom","SegmentStart","SegmentEnd","Segment"]
        annotation_df.set_index("Segment")
        cnv_by_clone.set_index("Segment")
        copy_number_data = annotation_df.join(cnv_by_clone, lsuffix="annot", rsuffix="cnv")
        copy_number_data = copy_number_data[["GeneSymbol"] + list(index.keys())]
        copy_number_data.drop_duplicates(inplace=True)
        gene_subset = set(copy_number_data["GeneSymbol"]).intersection(set(gene_subset))
        gene_subset = list(gene_subset)
        copy_number_data = copy_number_data[copy_number_data["GeneSymbol"].isin(gene_subset)]
        copy_number_data = copy_number_data.dropna()
        copy_number_data.to_csv("./tables/copy_number_data.tsv", sep="\t", index=False)
        copy_number_data = copy_number_data.set_index("GeneSymbol")
        return copy_number_data

    @staticmethod
    def write_input(gene_expression_data, copy_number_data, clone_assignments, genes):
        assembled_cnv = CloneAlign.assemble_copy_number_data(copy_number_data, clone_assignments, genes)
        genes = list(assembled_cnv.index)
        print(len(genes))
        gene_expression_data = gene_expression_data[:,genes].X.todense()
        valid_rows = []
        for row in gene_expression_data:
            row = numpy.asarray(row).tolist()[0]
            if len(row) != row.count(0.0):
                valid_rows.append(row)
        gene_expression_data = numpy.array(valid_rows)
        print(assembled_cnv.shape)
        print(gene_expression_data.shape)
        robjects.r.assign("gene_expression_matrix", gene_expression_data)
        robjects.r("saveRDS(gene_expression_matrix,file='rdata/gene_expression_data.rdata')")
        robjects.r.assign("copy_number_data", assembled_cnv)
        robjects.r("saveRDS(copy_number_data,file='rdata/copy_number_data.rdata')")

    @staticmethod
    def run_command_line(gene_expression_matrix, copy_number_data, clone_assignments, genes):
        CloneAlign.write_input(gene_expression_matrix, copy_number_data, clone_assignments, genes)
        CloneAlign.command()
        subprocess.call(["Rscript","run_clonealign.R"])

    @staticmethod
    def command():
        script = open("run_clonealign.R","w")
        script.write("""
        library(clonealign)
        reticulate::use_python("/home/nceglia/.virtualenvs/r-tensorflow/bin/python", required=TRUE)
        gene_expression_data <- readRDS("./rdata/gene_expression_data.rdata")
        gene_expression_data[is.na(gene_expression_data)] <- 0
        copy_number_data <- readRDS("./rdata/copy_number_data.rdata")
        cal <- clonealign(gene_expression_data, copy_number_data)
        saveRDS(cal,"./rdata/clone_align_fit.rdata")
        """)
        script.close()
