import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import pickle

from interface.tenxanalysis import TenxAnalysis
from software.cellassign import CellAssign
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

from utils.config import Configuration

def Run(sampleid, before, finished):
    open(finished,"w").write("Completed")


def RunClustering(sampleid, workflow):
    workflow.transform (
        name = "clustering",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("qc.complete"),
            pypeliner.managed.OutputFile("clustering.complete")
        )
    )
    return workflow
