import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
from multipledispatch import dispatch
from functools import partial
import os
import pandas
import numpy
import argparse
import warnings

warnings.filterwarnings("ignore")
pandas2ri.activate()

SingleCellExperimentInterface = importr('SingleCellExperiment')
SummarizedExperimentInterface = importr('SummarizedExperiment')
BiocGenericsInterface         = importr('BiocGenerics')


"""
SingleCellExperiment

Specific Fields:
    int_elementMetadata
    int_colData
    int_metadata
    reducedDims

Inherited from SummarizedExperiment:
    assays
    colData
    rowData

Inherited from BiocGenerics
    sizeFactors
"""

class SingleCellExperiment(RS4):

    @classmethod
    def fromRData(sce_class, rdata):
        rs4_object = robjects.r[r.load(rdata)[0]]
        return sce_class.fromRS4(rs4_object)

    @classmethod
    def fromRS4(sce_class, rs4_object):
        sce = sce_class(rs4_object)
        sce.rowData = SummarizedExperimentInterface.rowData(sce)
        sce.colData = SummarizedExperimentInterface.colData(sce)
        sce.assays = SummarizedExperimentInterface.assays(sce)
        sce.reducedDims = SingleCellExperimentInterface.reducedDims(sce)
        sce.sizeFactors = BiocGenericsInterface.sizeFactors(sce)
        return sce

    @staticmethod
    def unpack(rs4_object):
        unpacked_object = dict()
        if type(rs4_object) == rinterface.RNULLType:
            return unpacked_object
        for slot in rs4_object.slotnames():
            if type(rs4_object.slots[slot]) == rinterface.RNULLType: continue
            robj = pandas2ri.ri2py(rs4_object.slots[slot])
            if type(robj) == robjects.vectors.ListVector:
                list_array = list()
                for x in pandas2ri.ri2py(robj):
                    list_array.append(pandas2ri.ri2py(x))
                unpacked_object[slot] = list_array
            else:
                unpacked_object[slot] = list(robj)
        return unpacked_object

    @staticmethod
    def get_axis(axis_data):
        assert "nrows" in axis_data, "Axis requires nrows"
        assert "rownames" in axis_data, "Axis requires rownames"
        assert "listData" in axis_data, "Axis requires listData"
        nrows = axis_data["nrows"]
        rownames = axis_data["rownames"]
        categorical = axis_data["listData"]
        return rownames, nrows[0], categorical

    @property
    def assayNames(self):
        return tuple(self._assays.keys())

    def assay(self, assay):
        row_names, nrows, _ = self.rowData
        col_names, ncols, _ = self.colData
        assert (nrows, ncols) == self.assays[assay].shape, "Invalid Data Shape"
        row_names_index = pandas.Index(row_names, name="rows")
        col_names_index = pandas.Index(col_names, name="columns")
        dataframe = pandas.DataFrame(data=self.assays[assay],
                                     index=row_names_index,
                                     columns=col_names_index)
        return dataframe

    @property
    def rowData(self):
        return SingleCellExperiment.get_axis(self._rowData)

    @rowData.setter
    def rowData(self, rs4_rowData):
        self._rowData = SingleCellExperiment.unpack(rs4_rowData)

    @property
    def colData(self):
        return SingleCellExperiment.get_axis(self._colData)

    @colData.setter
    def colData(self, rs4_colData):
        self._colData = SingleCellExperiment.unpack(rs4_colData)

    @property
    def assays(self):
        return self._assays

    @assays.setter
    def assays(self, rs4_assays):
        list_vector = pandas2ri.ri2py(rs4_assays.slots["listData"])
        assays = dict()
        for assay, label in zip(list_vector, list_vector.names):
            if type(assay) == robjects.methods.RS4:
                assays[label] = dict()
                for slot in assay.slotnames():
                    assays[label][slot] = pandas2ri.ri2py(assay.slots[slot])
            elif type(assay) == robjects.vectors.Matrix:
                assays[label] = pandas2ri.ri2py(assay)
        self._assays = assays

if __name__ == '__main__':
    sce = SingleCellExperiment.fromRData("~/data/example_sce.RData")
    for assay in sce.assayNames:
        print(assay,sce.assay(assay))

    # matrices = sce.assayList
    # print(sce.rowData)
    # print(sce.colData)
    # print(sce.assays)
    # print(sce.reducedDims)
    # print(sce.sizeFactors)
