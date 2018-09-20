import argparse
import glob
import inspect
import os
import numpy
import pickle
from interface.singlecellexperiment import SingleCellExperiment
from interface.dropletutils import DropletUtils

def exportRData(rdata, directory, delim="\t"):
    if not os.path.exists(directory):
        os.makedirs(directory)
    sce = SingleCellExperiment.fromRData(rdata)
    output = open(os.path.join(directory,"meta.txt"),"w")
    output.write("sizeFactors: " + str(sce.sizeFactors) + "\n")
    output.write("reducedDims: " + str(sce.reducedDims) + "\n")
    for assay in sce.assayNames:
        filename = os.path.join(directory,"{}.csv".format(assay))
        print(filename)
        dataframe = sce.assay(assay)
        dataframe.to_csv(filename,sep=delim)

def exportRS4(rs4_object, directory, delim="\t"):
    if not os.path.exists(directory):
        os.makedirs(directory)
    output = open(os.path.join(directory,"meta.txt"),"w")
    output.write("sizeFactors: " + str(sce.sizeFactors) + "\n")
    output.write("reducedDims: " + str(sce.reducedDims) + "\n")
    for assay in rs4_object.assayNames:
        filename = os.path.join(directory,"{}.csv".format(assay))
        print(filename)
        dataframe = rs4_object.assay(assay)
        dataframe.to_csv(filename,sep=delim)


def inportRS4(rs4_directory, delim="\t"):
    raise AssertionError("Not Implemented")

if __name__ == '__main__':
    pass
