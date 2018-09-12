import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
from multipledispatch import dispatch
from functools import partial
import os
import pandas
import argparse
import warnings

warnings.filterwarnings("ignore")
pandas2ri.activate()

SingleCellExperimentInterface = importr('SingleCellExperiment')
SummarizedExperimentInterface = importr('SummarizedExperiment')
ScranInterface                = importr('scran')
BioBaseInterface              = importr('Biobase')
DESeq2Interface               = importr('DESeq2')
BiocGenerics                  = importr('BiocGenerics')


"""
SCE Specific Fields:
    int_elementMetadata
    int_colData
    int_metadata
    reducedDims


"""

class SingleCellExperiment(RS4):

    @classmethod
    def fromRData(sce_class, rdata):
        data = robjects.r[r.load(rdata)[0]]
        return sce_class(data)

    @classmethod
    def fromRS4(sce_class, rs4_object):
        return sce_class(rs4_object)

    def rowData(self):
        return SummarizedExperimentInterface.rowData(self)

    def colData(self):
        return SummarizedExperimentInterface.colData(self)

    def rowNames(self):
        return SummarizedExperimentInterface.rownames(self)

    def assayData(self):
        res = SummarizedExperimentInterface.assays(self)
        list_vector = pandas2ri.ri2py(res.slots["listData"])
        assays = dict()
        for assay, label in zip(list_vector, list_vector.names):
            if type(assay) == robjects.methods.RS4:
                assays[label] = dict()
                for slot in assay.slotnames():
                    assays[label][slot] = pandas2ri.ri2py(assay.slots[slot])
            elif type(assay) == robjects.vectors.Matrix:
                assays[label] = pandas2ri.ri2py(assay)
        return assays
        matrix = list()
        for x in list_vector.__dict___.keys():
            print(x)
        print(list_vector.names)
        print(dict(list_vector))
        json = dict(zip(list(list_vector.names), matrix ))
        dataframe = pandas.DataFrame.from_dict(json)
        return dataframe

    def reducedDims(self):
        return SingleCellExperimentInterface.reducedDims(self)

    def sizeFactors(self):
        return BiocGenerics.sizeFactors(self)


if __name__ == '__main__':
    sce = SingleCellExperiment.fromRData("~/data/example_sce.RData")
    print(sce.rowData())
    print(sce.colData())
    print(sce.assayData())
    print(sce.reducedDims())
    print(sce.sizeFactors())
