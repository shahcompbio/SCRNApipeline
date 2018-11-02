import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr

import pandas
import numpy

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
