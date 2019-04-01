import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

"""cellranger count --id=TENX_014_SA604X6XB01979_001 --mempercore=8 --transcriptome=/work/shah/ceglian/refdata/mouse --fastqs=/juno/work/shah/ceglian/temperature/TENX_014_SA604X6XB01979_001 --maxjobs=200 --jobmode=lsf"""

def Run(sampleid, mito):
    fastq = FastQDirectory(fastq_directory, self.prefix, self.output, datapath=config.datapath)
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis,sampleid)
    qc.run()
    qc.upload()

def RunCellranger(fastqs, workflow):
    self.workflow.transform (
        name = "{}_cellranger_counts".format(self.prefix),
        func = CellRanger.count,
        args = (
            self.fastqs,
            pypeliner.managed.OutputFile()
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    mito = sys.argv[2]
    Run(sampleid, mito)
