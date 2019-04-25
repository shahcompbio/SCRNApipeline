import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from software.cellranger import CellRanger
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage

from utils.config import Configuration, write_config

config = Configuration()

def Run(sampleid, finished):
    CellRanger.count([sampleid])
    open(finished,"w").write("Completed")

def upload_tenx(sampleid, before, finished):
    print("Calling upload.")
    # tenx = TenxAnalysis("./{}/outs/".format(sampleid))
    # tenx.finalize()
    # tenxds = TenxDataStorage(sampleid)
    # tenxds.upload_cellranger(tenx)
    open(finished,"w").write("Completed")

def RunCellranger(sampleid, workflow):
    workflow.transform (
        name = "cellranger_counts",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.OutputFile("cellranger.complete")
        )
    )
    workflow.transform (
        name = "upload_tenx",
        func = upload_tenx,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellranger.complete"),
            pypeliner.managed.OutputFile("upload.complete")
        )
    )
    return workflow
