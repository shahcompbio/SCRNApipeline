import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from interface.tenxanalysis import TenxAnalysis
from software.cellassign import CellAssign
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

from utils.config import Configuration

config = Configuration()

def Run(sampleid, before, finished):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    analysis_path = tenx.tenx_path
    tenx_analysis = TenxAnalysis(analysis_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    qc.build()
    fit = os.path.join(analysis_path, "{}_cellasign.rdata".format(sampleid))
    CellAssign.run(qc.sce, config.rho_matrix, fit)
    open(finished,"w").write("Completed")


def RunCellAssign(sampleid, workflow):
    workflow.transform (
        name = "cellassign",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("qc.complete"),
            pypeliner.managed.OutputFile("cellassign.complete")
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid, "qc.complete", "cellassign.complete")
