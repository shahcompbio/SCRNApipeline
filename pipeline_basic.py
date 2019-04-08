import argparse
import glob
import os
import yaml
import sys

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from utils.config import Configuration, write_config

from workflows.run_cellranger import RunCellranger
from workflows.run_qc import RunQC

def create_workflow():

    workflow = pypeliner.workflow.Workflow()

    _fastqs = args.get("sampleid")
    prefix = args.get("sampleid")
    reference = args.get("reference")
    output = args.get("output")
    write_config(prefix,output,reference)

    print("Setting up fastqs...")
    fastqs = []
    fastqs.append(FastQDirectory(fastq, prefix, output, "/data"))

    workflow = RunCellranger(fastqs, workflow)

    return workflow

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--sampleid', type=str, nargs="+", help='Sample ID linked to fastqs in scrnadata.')
    argparser.add_argument('--output', type=str, help="Base directory for output")
    argparser.add_argument("--reference", type=str, help="Reference build.")

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
