import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr

import pandas

# SingleCellExperimentInterface = importr('')

class GeneMarkerMatrix(object):

    def __init__(self):
        pass

    def load_example(self):
        pandas2ri.activate()
        matrix = robjects.r[r.load("./tests/gene_marker_matrix.rdata")[0]]
        dataframe = []
        gene_names = ["Gene{}".format(i) for i in range(10)]
        for gene, row in zip(gene_names, matrix):
            dataframe.append((gene, list(map(int,list(row)))))
        dataframe = pandas.DataFrame.from_items(dataframe, orient="index", columns=["Cell1","Cell2"])
        dataframe = pandas2ri.py2ri(dataframe)
        return robjects.r.matrix(dataframe)
