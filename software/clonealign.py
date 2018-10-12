import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.rinterface import RRuntimeError

import string
import numpy
import subprocess
import os
import warnings
import glob
import sys

from interface.singlecellexperiment import SingleCellExperiment
from interface.clonealignfit import CloneAlignFit

warnings.filterwarnings("ignore")

os.environ["RETICULATE_PYTHON"] = sys.executable

CloneAlignInterface = importr("clonealign")

class CloneAlign(object):


    @staticmethod
    def run(sce_experiment, cnv_data):
        assert "counts" in sce_experiment.assayNames, "No counts in assays"
        matrix = numpy.transpose(sce_experiment.assays["counts"].toarray())
        cal = CloneAlignInterface.clonealign(matrix, cnv_data)
        cal = dict(zip(cal.names, map(list,list(cal))))
        return cal
