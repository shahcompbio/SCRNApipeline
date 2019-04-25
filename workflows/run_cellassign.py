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
from interface.genemarkermatrix import GeneMarkerMatrix
from utils.plotting import celltypes, tsne_by_cell_type, umap_by_cell_type

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
    CellAssign.run(qc.sce, config.rho_matrix, ".cache/{}/celltypes.rdata".format(sampleid))
    open(finished,"w").write("Completed")

def Analysis(sampleid, before, finished):
    tenx = TenxDataStorage(sampleid, version="v3")
    tenx.download()
    analysis_path = tenx.tenx_path
    tenx_analysis = TenxAnalysis(analysis_path)
    tenx_analysis.load()
    tenx_analysis.extract()
    qc = QualityControl(tenx_analysis, sampleid)
    cellassign_analysis = ".cache/{}/cellassignanalysis/".format(sampleid)
    if not os.path.exists(cellassign_analysis):
        os.makedirs(cellassign_analysis)
    pyfit = os.path.join(".cache/{}/cell_types.pkl".format(sampleid))
    assert os.path.exists(pyfit), "No Pyfit Found."
    pyfit = pickle.load(open(pyfit,"rb"))
    marker_list = GeneMarkerMatrix.read_yaml(config.rho_matrix)
    cell_types = marker_list.celltypes()
    if "B cell" not in cell_types: cell_types.append("B cell")
    celltypes(pyfit, sampleid, cellassign_analysis, known_types=cell_types)
    tsne_by_cell_type(qc.sce, pyfit, sampleid, cellassign_analysis, known_types=cell_types)
    umap_by_cell_type(qc.sce, pyfit, sampleid, cellassign_analysis, known_types=cell_types)
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
    workflow.transform (
        name = "cellassignanalysis",
        func = Analysis,
        args = (
            sampleid,
            pypeliner.managed.InputFile("cellassign.complete"),
            pypeliner.managed.OutputFile("cellassignanalysis.complete")
        )
    )
    return workflow

if __name__ == '__main__':
    sampleid = sys.argv[1]
    Run(sampleid, "qc.complete", "cellassign.complete")
