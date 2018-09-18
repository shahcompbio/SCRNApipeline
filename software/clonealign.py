import urllib.request
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from urllib.error import HTTPError
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr
import string
import numpy
import pandas

MatrixInterface               = importr('Matrix')

class CloneAlign(object):

    def __init__(self):
        decode = lambda x: x.decode("utf-8")
        try:
            github_url = urllib.request.urlopen('https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/R/clonealign.R')
            code = '\n'.join(map(decode, github_url.readlines()))
            self.clonealign = SignatureTranslatedAnonymousPackage(code, "clonealign")
        except HTTPError:
            pass

    @staticmethod
    def loadFitMatrix(rdata_path):
        rs4_obj = robjects.r[r.load(rdata_path)[0]]
        for field in rs4_obj:
            if type(field) == rinterface.RNULLType: continue
            if type(field) == robjects.vectors.ListVector:
                for col in field:
                    print("Coltype",type(col))
        return None

    def dumpFitMatrix(rdata_path):
        pass

    def run(self, sce_experiment, copy_number_data):
        #copy_number_matrix = CloneAlign.loadMatrix(copy_number_data)
        _, nrows, data = sce_experiment.rowData
        _, ncols, _ = sce_experiment.colData
        print(nrows, "nrows")
        print(ncols, "ncols")
        df = dict()
        for col, row in zip(list(string.ascii_uppercase), data):
            df[col] = row
        copy_number_matrix = pandas2ri.py2ri(pandas.DataFrame.from_dict(df))
        return self.clonealign.clonealign(sce_experiment, copy_number_matrix)
