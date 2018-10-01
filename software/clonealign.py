import string
import numpy
import subprocess
import os
import warnings
import glob

from interface.singlecellexperiment import SingleCellExperiment
from interface.clonealignfit import CloneAlignFit

warnings.filterwarnings("ignore")

cwd = os.path.dirname(os.path.realpath(__file__))
tmp_assay = "software/scripts/tmp_assay.tsv"
tmp_alignments = "scripts/tmp_clone_alignments.tsv"
tmp_params = "scripts/tmp_params.tsv"

class CloneAlign(object):

    @staticmethod
    def clean_up():
        tmps = glob.glob(os.path.join(cwd, "scripts/tmp_*"))
        print(tmps)

    @staticmethod
    def call():
        script = os.path.join(cwd, "scripts/run_clonealign.R")
        return subprocess.call(["Rscript", script])

    @staticmethod
    def run(sce_experiment):
        assert "counts" in sce_experiment.assayNames, "No counts in assays"
        counts = sce_experiment.assays["counts"]
        numpy.savetxt(tmp_assay, counts.todense(), delimiter="\t")
        ret = CloneAlign.call()
        assert ret == 0
        assert os.path.exists(os.path.join(cwd, tmp_alignments))
        fit = CloneAlignFit(clones=os.path.join(cwd, tmp_alignments),
                            genes=sce_experiment.rowData[0])
        CloneAlign.clean_up()
        return fit
