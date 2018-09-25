import urllib.request
from interface.singlecellexperiment import SingleCellExperiment
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
import gc

class CloneAlign(object):

    @staticmethod
    def load_methods():
        # decode = lambda x: x.decode("utf-8")
        # github_url = urllib.request.urlopen('https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/R/clonealign.R')
        # code = '\n'.join(map(decode, github_url.readlines()))
        # github_url = urllib.request.urlopen('https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/R/inference-tflow.R')
        # code += '\n'.join(map(decode, github_url.readlines()))
        # stop_code = "stop(msg)"
        # if stop_code in code:
        #     print("Found stop code")
        #     code = code.replace(stop_code,"")
        # code = "library('tensorflow')\n" + code
        # code = "library('glue')\n" + code
        # code = "library('reticulate')\n" + code
        code = open('/home/ceglian/codebase/clonealign/R/clonealign.R',"r").read()
        return SignatureTranslatedAnonymousPackage(code, "clonealign")

    @staticmethod
    def loadFitMatrix(rdata_path):
        rs4_obj = robjects.r[r.load(rdata_path)[0]]
        for field in rs4_obj:
            if type(field) == rinterface.RNULLType: continue
            if type(field) == robjects.vectors.ListVector:
                for col in field:
                    print("Coltype",type(col))
        return None

    @staticmethod
    def dumpFitMatrix(rdata_path):
        pass

    @staticmethod
    def run(sce_experiment):
        # ReticulateInterface = importr('reticulate')
        # TensorFlowInterface = importr('tensorflow')
        # #TensorFlowInterface.install_tensorflow()
        # CloneAlignInterface = importr('clonealign')
        clonealign = CloneAlign.load_methods()
        rownames, nrows, data = sce_experiment.rowData
        _, ncols, _ = sce_experiment.colData
        df = dict()
        columns = []
        for col, row in zip(list(string.ascii_uppercase), data):
            df[col] = row
            columns.append(col)
        df["names"] = rownames
        df = pandas.DataFrame.from_dict(df)
        df.set_index(["names"])
        copy_number_matrix = pandas2ri.py2ri(df)
        return clonealign.clonealign(sce_experiment, copy_number_matrix, max_iter_em=5)
