import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr

import os
import pandas
import numpy
import argparse
import warnings
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")
pandas2ri.activate()

SingleCellExperimentInterface = importr('SingleCellExperiment')
SummarizedExperimentInterface = importr('SummarizedExperiment')
BiocGenericsInterface         = importr('BiocGenerics')
MatrixInterface               = importr('Matrix')
BiobaseInterface              = importr('Biobase')

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

    def __eq__(self, other):
        # TODO better assertion of equal objects
        return self.assayNames == other.assayNames

    @classmethod
    def fromRData(sce_class, rdata):
        # try:
        #     rs4_object = robjects.r[r.load(rdata)[0]]
        # except Exception as e:
        rs4_object = r.readRDS(rdata)
        sce = sce_class.fromRS4(rs4_object)
        sce.rs4 = rs4_object
        return sce

    def asSummarizedExperiment(self):
        data = robjects.r["SummarizedExperiment"](self.rs4)
        data.slots["assays"] = SummarizedExperimentInterface.assays(self.rs4)
        return data

    def asExpressionSet(self):
        sme = self.asSummarizedExperiment()
        data = robjects.r["ExpressionSet"](assayData=sme.slots["assays"])
        return data

    def save(self,filename):
        robjects.r.assign("sce",self.rs4)
        robjects.r("saveRDS(sce, file='{}')".format(filename))

    @classmethod
    def fromRS4(sce_class, rs4_object):
        sce = sce_class(rs4_object)
        sce.rs4 = rs4_object
        sce.rowData = SummarizedExperimentInterface.rowData(sce)
        sce.colData = SummarizedExperimentInterface.colData(sce)
        sce.assays = SummarizedExperimentInterface.assays(sce)
        sce.reducedDims = SingleCellExperimentInterface.reducedDims(sce)
        sce.sizeFactors = BiocGenericsInterface.sizeFactors(sce)
        sce.rownames = robjects.r["rownames"](sce)
        sce.colnames = robjects.r["colnames"](sce)
        return sce

    @classmethod
    def toSummarizedExperiment(sce_class, rs4_object):
        return SummarizedExperimentInterface.SummarizedExperiment()

    @staticmethod
    def unpack(rs4_object):
        unpacked_object = dict()
        if type(rs4_object) == rinterface.RNULLType:
            return unpacked_object
        for slot in rs4_object.slotnames():
            value = rs4_object.slots[slot]
            if type(value) == rinterface.RNULLType:
                continue
            elif type(value) == robjects.vectors.ListVector:
                if value.names:
                    value = dict(zip(value.names, map(list,list(value))))
                    for column, data in value.items():
                        unpacked_object[column] = data
            else:
                unpacked_object[slot] = list(value)[0]
        return unpacked_object

    @property
    def assayNames(self):
        return tuple(self._assays.keys())

    def assay(self, assay, row_index=0, col_index=0):
        assay = self.assays[assay]
        row_names, nrows, row_data = self.rowData
        col_names, ncols, col_data = self.colData
        assert (nrows, ncols) == assay.shape
        row_index = pandas.MultiIndex.from_tuples(list(zip(*row_data)))
        col_index = pandas.MultiIndex.from_tuples(list(zip(*col_data)))
        assay_df = pandas.DataFrame(data=assay,index=row_index,columns=col_index)
        return assay_df

    @property
    def rowData(self):
        return self._rowData

    @rowData.setter
    def rowData(self, rs4_rowData):
        self._rowData = SingleCellExperiment.unpack(rs4_rowData)

    @property
    def colData(self):
        return self._colData

    @colData.setter
    def colData(self, rs4_colData):
        self._colData = SingleCellExperiment.unpack(rs4_colData)

    @property
    def assays(self):
        return self._assays

    @staticmethod
    def DCGtoCSR(data, row_ind, col_pointers, nrows):
        data = list(data)
        begin_pointer = col_pointers[0]
        col_ind = numpy.zeros(len(row_ind))
        for column, pointer in enumerate(col_pointers[1:]):
            for row in range(begin_pointer, pointer):
                col_ind[row] = column
            begin_pointer = pointer
        return csr_matrix((data,(row_ind,col_ind)),shape=(nrows,len(col_pointers)-1))

    @staticmethod
    def CSRtoDCG(sparse_matrix):
        data = robjects.DataFrame(sparse_matrix.toarray().flatten())
        nrows, ncols = sparse_matrix.shape
        return MatrixInterface.Matrix(data, nrow=nrows, ncol=ncols, sparse=True)

    @assays.setter
    def assays(self, rs4_assays):
        list_vector = pandas2ri.ri2py(rs4_assays.slots["listData"])
        self._assays = dict()
        for assay, label in zip(list_vector, list_vector.names):
            if type(assay) == robjects.methods.RS4:
                non_zero_elements = assay.slots["x"]
                row_numbers =pandas2ri.ri2py(assay.slots["i"])
                column_pointers = pandas2ri.ri2py(assay.slots["p"])
                nrows = len(list(pandas2ri.ri2py(assay.slots["Dimnames"]))[0])
                self._assays[label] = SingleCellExperiment.DCGtoCSR(non_zero_elements, row_numbers, column_pointers, nrows)
            elif type(assay) == robjects.vectors.Matrix:
                self._assays[label] = csr_matrix(pandas2ri.ri2py(assay))



if __name__ == '__main__':
    pass
