import argparse
import glob
import os
import yaml
import sys

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from utils.config import Configuration

from workflows.run_cellranger import RunCellranger
from workflows.run_qc import RunQC
from workflows.run_report import RunReport
from workflows.run_cellassign import RunCellAssign

config = Configuration()

def create_workflow():
    print("Creating workflow.")
    workflow = pypeliner.workflow.Workflow()

    prefix = args['sampleid']

    workflow = RunCellranger(prefix, workflow)
    workflow = RunQC(prefix, workflow)
    workflow = RunReport(prefix, workflow)

    return workflow

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
