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

def BuildInput(sampleid, before, finished):
    open(finished,"w").write("Completed")

def Run(sampleid, before, finished):
    open(finished,"w").write("Completed")

def Analysis(sampleid, before, finished):
    open(finished,"w").write("Completed")


def RunScvis(sampleid, workflow):
    workflow.transform (
        name = "scvis_input",
        func = BuildInput,
        args = (
            sampleid,
            pypeliner.managed.InputFile("qc.complete"),
            pypeliner.managed.OutputFile("scvis_input.complete")
        )
    )
    workflow.transform (
        name = "run_scvis",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("scvis_input.complete"),
            pypeliner.managed.OutputFile("scvis.complete")
        )
    )
    workflow.transform (
        name = "scvis_analysis",
        func = Analysis,
        args = (
            sampleid,
            pypeliner.managed.InputFile("scvis.complete"),
            pypeliner.managed.OutputFile("scvis_analysis.complete")
        )
    )

    return workflow
