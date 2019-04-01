import argparse
import glob
import os
import yaml
import sys

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from software.cellranger import CellRanger
from software.cellassign import CellAssign
from software.clonealign import CloneAlign
from software.scviz import SCViz

from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from interface.genemarkermatrix import GeneMarkerMatrix, generate_json
from interface.singlecellexperiment import SingleCellExperiment

from utils.config import Configuration
from utils.cloud import TenxDataStorage, VirtualMachine

from workflows.run_cellranger import RunCellranger
from workflows.run_qc import RunQC



config = Configuration()

def create_workflow():

    workflow = pypeliner.workflow.Workflow()


    RunCellranger()
    RunQC()




if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--fastqs', type=str, nargs="+", help='CellRanger Structured FastQ Output Directory')
    argparser.add_argument('--out', type=str, help="Base directory for output")
    argparser.add_argument("--prefix", type=str, help="Analysis prefix")
    argparser.add_argument("--yaml", type=str, help="Configuration settings for pipeline.")

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
