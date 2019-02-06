import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
import scanpy.api as sc
import collections
import operator
import pandas
import numpy
import os

def generate_json(analysis, sce, organ):
    markers = os.path.join(os.path.dirname(os.path.realpath(__file__)),"../utils/marker_genes.tsv")
    rows = open(markers,"r").read().splitlines()
    header = rows.pop(0).split("\t")
    matrix = collections.defaultdict(list)
    for row in rows:
        row = dict(zip(header,row.split("\t")))
        if "Hs" in row["species"] and row["organ"].lower() == organ.lower():
            if float(row["ubiquitousness index"]):
                matrix[row["cell type"]].append(row["official gene symbol"])
                aliases = row["nicknames"].split("|")
                for alias in aliases:
                    matrix[row["cell type"]].append(alias)
    markers = analysis.markers(sce)
    adata = analysis.create_scanpy_adata(sce, high_var = True)
    final_set = list(set(adata.var.index).intersection(set(markers)))
    adata = adata[:,final_set]
    rho = collections.defaultdict(list)
    for cell_type in matrix.keys():
        count = 0
        for gene in matrix[cell_type]:
            if gene in adata.var_names:
                vals = adata[:, gene].X
                dropout = list(vals).count(0)
                if dropout == adata.shape[0]:
                    continue
                rho[cell_type].append(gene)
    return rho

class GeneMarkerMatrix(object):

    def __init__(self, marker_list):
        self.marker_list = marker_list
        self.genes = list(set([gene for genelist in self.marker_list.values() for gene in genelist]))
        self.cells = list(self.marker_list.keys())

    def matrix(self, subset=None, include_other=True):
        self.indicator = []
        if subset != None:
            self.genes = list(set([gene for genelist in self.marker_list.values() for gene in genelist]).intersection(subset))
        else:
            self.genes = list(set([gene for genelist in self.marker_list.values() for gene in genelist]))
        self.cells = list(self.marker_list.keys())
        if include_other:
            self.cells.append("Other")
        for cell in self.cells:
            marker_genes = []
            for gene in self.genes:
                if cell == "Other":
                    marker_genes.append(0)
                elif gene in self.marker_list[cell]:
                    marker_genes.append(1)
                else:
                    marker_genes.append(0)
            self.indicator.append(marker_genes)
        return numpy.transpose(numpy.matrix(self.indicator))

    def celltypes(self):
        return self.cells

    def matrix_from_rdata(self, rdata):
        pandas2ri.activate()
        matrix = robjects.r[r.load("./tests/gene_marker_matrix.rdata")[0]]
        return matrix
