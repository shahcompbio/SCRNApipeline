import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import numpy
import os
import sys

os.environ["RETICULATE_PYTHON"] = sys.executable

class CloneAlign(object):

    @staticmethod
    def run(sce_experiment, cnv_data, filename):
        CloneAlignInterface = importr("clonealign")
        sce_experiment = SingleCellExperiment.fromRData(sce_experiment)
        assert "counts" in sce_experiment.assayNames, "No counts in assays"
        matrix = numpy.transpose(sce_experiment.assays["counts"].toarray())
        cal = CloneAlignInterface.clonealign(matrix, cnv_data)
        robjects.r.assign("clone_align_fit", clone_align_fit)
        robjects.r("saveRDS(clone_align_fit, file='{}')".format(filename))
