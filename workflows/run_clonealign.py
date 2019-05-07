import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os
import pickle

from interface.tenxanalysis import TenxAnalysis
from software.clonealign import CloneAlign
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

from utils.config import Configuration

def BuildInput(dlpid, before, finished):
    open(finished,"w").write("Completed")

def Run(sampleid, before, finished):
    open(finished,"w").write("Completed")

def Analysis(sampleid, before, finished):
    open(finished,"w").write("Completed")

def RunCellAssign(sampleid, workflow):
    workflow.transform (
        name = "builddlpinput",
        func = BuildInput,
        args = (
            sampleid,
            pypeliner.managed.InputFile("qc.complete"),
            pypeliner.managed.OutputFile("builddlp.complete")
        )
    )
    workflow.transform (
        name = "clonealign",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("builddlp.complete"),
            pypeliner.managed.OutputFile("clonealign.complete")
        )
    )
    workflow.transform (
        name = "clonealignanalysis",
        func = Analysis,
        args = (
            sampleid,
            pypeliner.managed.InputFile("clonealign.complete"),
            pypeliner.managed.OutputFile("clonealignanalysis.complete")
        )
    )
    return workflow
