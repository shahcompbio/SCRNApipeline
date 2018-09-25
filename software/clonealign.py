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
import subprocess
import os

import warnings
warnings.filterwarnings("ignore")


tmp_matrix = "./software/scripts/tmp_matrix.tsv"

class CloneAlign(object):

    # @staticmethod
    # def load_methods():
    #     # decode = lambda x: x.decode("utf-8")
    #     # github_url = urllib.request.urlopen('https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/R/clonealign.R')
    #     # code = '\n'.join(map(decode, github_url.readlines()))
    #     # github_url = urllib.request.urlopen('https://raw.githubusercontent.com/kieranrcampbell/clonealign/master/R/inference-tflow.R')
    #     # code += '\n'.join(map(decode, github_url.readlines()))
    #     # stop_code = "stop(msg)"
    #     # if stop_code in code:
    #     #     print("Found stop code")
    #     #     code = code.replace(stop_code,"")
    #     # code = "library('tensorflow')\n" + code
    #     # code = "library('glue')\n" + code
    #     # code = "library('reticulate')\n" + code
    #     code = open('/home/ceglian/codebase/clonealign/R/clonealign.R',"r").read()
    #     return SignatureTranslatedAnonymousPackage(code, "clonealign")

    @staticmethod
    def clean_up():
        tmps = glob.glob("./scripts/tmp_*")

    @staticmethod
    def call():
        cwd = os.path.dirname(os.path.realpath(__file__))
        script = os.path.join(cwd, "scripts/run_clonealign.R")
        output = os.path.join(cwd, "scripts/tmp_matrix.txt")
        subprocess.call(["Rscript", script, "-f", output])

    @staticmethod
    def run(sce_experiment):
        assert "counts" in sce_experiment.assayNames, "No counts in assays"
        counts = sce_experiment.assays["counts"]
        numpy.savetxt(tmp_matrix, counts.todense(), delimiter="\t")
        ret = CloneAlign.call()
        print(ret)
        assert ret == 0
