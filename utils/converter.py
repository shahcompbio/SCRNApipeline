import rpy2.robjects as robjects
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
import os
import pandas
import argparse
import warnings

warnings.filterwarnings("ignore")

SingleCellExperimentInterface = importr('SingleCellExperiment')
pandas2ri.activate()

class SingleCellExperiment(RS4):

    def __init__(self, rdata_path):
        data = r.load(rdata_path)
        self.data = robjects.r[data[0]]
        for slot, value in self.data.slots.items():
            setattr(self, slot, value)

    def attributes(self):
        return list(self.data.slotnames())

    def matrix(self,index="Cell"):
        list_vector = pandas2ri.ri2py(self.colData.slots["listData"])
        columns = list_vector.names
        json = dict(zip(columns, map(list,list_vector)))
        dataframe = pandas.DataFrame.from_dict(json)
        dataframe.set_index(index, inplace=True)
        return dataframe

    def dump(self,filename,sep="\t"):
        sce.matrix().to_csv(filename,sep=sep)

    def info(self):
        info = dict()
        for key, value in self.colData.slots.items():
            if type(value) == robjects.vectors.IntVector:
                info[key] = list(value)
            elif type(value) == robjects.vectors.StrVector:
                info[key] = str(value)
            elif type(value) == robjects.vectors.ListVector:
                info[key] = list(value)
            else:
                continue
        return info


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Dump Single Cell Experiment RData to TSV')
    parser.add_argument('--rdata', type=str, dest='rdata', help='Input RData Single Cell Experiment')
    parser.add_argument('--output',type=str, dest='output',help="Output TSV")
    args = parser.parse_args()
    sce = SingleCellExperiment(args.rdata)
    sce.dump(args.output)
