import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr

import pandas
import numpy

# SingleCellExperimentInterface = importr('')

class GeneMarkerMatrix(object):

    def __init__(self, marker_list):
        self.indicator = []
        self.genes = list(set([gene for genelist in marker_list.values() for gene in genelist]))
        self.cells = marker_list.keys()
        for cell in self.cells:
            marker_genes = []
            for gene in self.genes:
                if gene in marker_list[cell]:
                    marker_genes.append(1)
                else:
                    marker_genes.append(0)
            self.indicator.append(marker_genes)
        self.indicator = numpy.matrix(self.indicator)

    def matrix(self):
        return numpy.transpose(self.indicator)

    def matrix_from_rdata(self, rdata):
        pandas2ri.activate()
        matrix = robjects.r[r.load("./tests/gene_marker_matrix.rdata")[0]]
        print("RHO",matrix.shape)
        return matrix
        # dataframe = []
        # for gene, row in zip(self.genes, matrix):
        #     dataframe.append((gene, list(map(int,list(row)))))
        # dataframe = pandas.DataFrame.from_items(dataframe, orient="index", columns=["Group1","Group2"])
        # dataframe = pandas2ri.py2ri(dataframe)
        # return robjects.r.matrix(dataframe)
        # print(ret)
        # return ret
