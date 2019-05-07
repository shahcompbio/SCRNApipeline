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

def RunClusteringDE(sample,before,finished):
    open(finished,"w").write("Completed")

def RunCellTypeDE(sample,before,finished):
    open(finished,"w").write("Completed")

def RunCloneDE(sample,before,finished):
    open(finished,"w").write("Completed")

def RunDifferentialAnalysis(sampleid, workflow):
    workflow.transform (
        name = "clustering_de",
        func = RunClusteringDE,
        args = (
            sampleid,
            pypeliner.managed.InputFile("clustering.complete"),
            pypeliner.managed.OutputFile("clustering_de.complete")
        )
    )
    workflow.transform (
        name = "celltype_de",
        func = RunCellTypeDE,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellassign.complete"),
            pypeliner.managed.OutputFile("celltype_de.complete")
        )
    )
    workflow.transform (
        name = "clone_de",
        func = RunCloneDE,
        args = (
            sampleid,
            pypeliner.managed.InputFile("clonealign.complete"),
            pypeliner.managed.OutputFile("clone_de.complete")
        )
    )

    return workflow
