import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys

from interface.tenxanalysis import TenxAnalysis
from software.cellassign import CellAssign
from utils.cloud import TenxDataStorage

from utils.config import Configuration

config = Configuration()

def Run(sampleid):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    sce = tenx_analysis.qcdrdata
    fit = os.path.join(tenx.tenx_path, "{}_cellasign.pkl".format(sampleid))
    rho_matrix = eval(open(config.rho_matrix,"r").read())
    CellAssign.run_em(tenx, sce, fit, sampleid, rho_matrix, None)


def RunCellAssign(sampleid):
    self.workflow.transform (
        name = "cellassign",
        func = CellAssign.run_em,
        args = (
            sampleid,
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid)
