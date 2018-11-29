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

def generate_json(analysis):
    rows = open("/home/nceglia/codebase/refdata/marker_genes.tsv","r").read().splitlines()
    header = rows.pop(0).split("\t")
    matrix = collections.defaultdict(list)
    for row in rows:
        row = dict(zip(header,row.split("\t")))
        if "Hs" in row["species"]:
            if float(row["ubiquitousness index"]) > 0.07:
                matrix[row["cell type"]].append(row["official gene symbol"])
                aliases = row["nicknames"].split("|")
                for alias in aliases:
                    matrix[row["cell type"]].append(alias)
    adata = analysis.create_scanpy_adata()
    adata.var_names_make_unique()
    sc.pp.recipe_zheng17(adata)
    rho = collections.defaultdict(list)
    types = dict()
    for cell_type in matrix.keys():
        count = 0
        for gene in matrix[cell_type]:
            if gene in adata.var_names:
                vals = adata[:, gene].X
                dropout = list(vals).count(0)
                #print(dropout, adata.shape)
                if dropout == adata.shape[0]:
                    print("dropped")
                    continue
                count += 1
        if count > 2:
            types[cell_type] = count
            rho[cell_type].append(gene)
    sorted_cells = sorted(types.items(), key=operator.itemgetter(1))
    for cell in sorted_cells:
        print(cell[0],cell[1])
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
