import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl
from utils.reporting import Results
from utils.export import exportMD, exportUpload

from utils.config import Configuration, write_config

config = Configuration()

def Run(sampleid, before, finished):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis,sampleid)
    plots = qc.plots
    cellassign = os.path.join(os.path.split(plots)[0],"cellassignanalysis")
    results = Results(config.jobpath)

    results.add_analysis(tenx.tenx_path)
    results.add_sce(qc.qcdsce)

    umi = os.path.join(plots,"umi.png")
    mito = os.path.join(plots,"mito.png")
    ribo = os.path.join(plots, "ribo.png")
    total_counts = os.path.join(plots, "total_counts.png")
    tfbc = os.path.join(plots, "total_features_by_counts.png")
    tcvfc = os.path.join(plots, "total_counts_v_features_by_counts.png")
    celltypes = os.path.join(cellassign, "cell_types.png")

    results.add_plot(umi,"UMI Distribution")
    results.add_plot(mito,"Mito Distribution")
    results.add_plot(ribo,"Ribo Distribution")
    results.add_plot(total_counts,"Total Counts Distribution")
    results.add_plot(tcvfc,"Total Counts")
    results.add_plot(tcvfc,"Total Features by Counts")
    results.add_plot(celltypes,"Cell Types")

    exportMD(results)
    exportUpload(results)
    open(finished,"w").write("Completed")


def RunReport(sampleid, workflow):
    workflow.transform (
        name = "report",
        func = Run,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellassignanalysis.complete"),
            pypeliner.managed.OutputFile("report.complete")
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid, None, "qc.complete")
