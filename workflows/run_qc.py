import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

def Run(sampleid, mito):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis,sampleid)
    qc.run()
    qc.upload()


def RunQC(tenx, mito, workflow):
    self.workflow.transform (
        name = "read10xcounts",
        func = Run,
        args = (
            tenx,
            pypeliner.managed.OutputFile(tenx.qcdsce)
        )
    )

if __name__ == '__main__':
    sampleid = sys.argv[1]
    mito = sys.argv[2]
    Run(sampleid, mito)
