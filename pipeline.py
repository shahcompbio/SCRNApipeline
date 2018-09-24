import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from software.cellranger import CellRanger
from software.tenx import TenX
from software.clonealign import CloneAlign
from software.cellassign import CellAssign

from interface.binarybasecall import BinaryBaseCall
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis

from utils.reporting import HTMLResults

import argparse

def create_workflow():

    bcl_directory = args["bcl_directory"]
    #fastq_directory = args["fastq_directory"]

    bcl_object = BinaryBaseCall(bcl_directory)
    #fastq_object = FastQDirectory(fastq_directory)
    #tenx_analysis = TenxAnalysis(fastq_object.out())

    workflow = pypeliner.workflow.Workflow()

    workflow.transform (
        name = "bcl_to_fastq",
        func = CellRanger.mkfastq,
        ret = pypeliner.managed.TempOutputObj("fastq_object"),
        args = (
            bcl_object,
        )
    )

    workflow.transform (
        name = "fastq_counts",
        func = CellRanger.count,
        ret = pypeliner.managed.TempOutputObj("tenx_analysis"),
        args = (
            pypeliner.managed.TempInputObj("fastq_object"),
        )
    )

    workflow.transform (
        name = "tenx_read10xcounts",
        func = TenX.read10xCounts,
        ret = pypeliner.managed.TempOutputObj("single_cell_experiment"),
        args = (
            pypeliner.managed.TempInputObj("tenx_analysis"),
        )
    )

    workflow.transform (
        name = "tenx_barcoderanks",
        func = TenX.barcodeRanks,
        args = (
            pypeliner.managed.TempInputObj("tenx_analysis"),
        )
    )

    workflow.transform (
        name = "tenx_emptydrops",
        func = TenX.emptyDrops,
        args = (
            pypeliner.managed.TempInputObj("tenx_analysis"),
        )
    )

    workflow.transform (
        name = "clonealign",
        func = CloneAlign.run,
        ret = pypeliner.managed.TempOutputObj("clone_align_fit")
        args = (
            pypeliner.managed.TempInputObj("single_cell_experiment")
        )
    )

    workflow.transform (
        name = "cellasign",
        func = CellAssign.run_em,
        ret = pypeliner.managed.TempOutputObj("cell_assignments"),
        args = (
            pypeliner.managed.TempInputObj("single_cell_experiment")
        )
    )

    workflow.transform (
        name = "html_output",
        func = HTMLResults.generate,
        args = (
            pypeliner.managed.TempInputObj("fastq_object")
            pypeliner.managed.TempInputObj("ten_analysis"),
            pypeliner.managed.TempInputObj("single_cell_experiment"),
            pypeliner.managed.TempInputObj("clone_align_fit"),
            pypeliner.managed.TempInputObj("cell_assignments"),
        )
    )

    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)
    argparser.add_argument('--bcl_directory', type=str, help='BaseCalls Illumina Directory')
    argparser.add_argument('--fastq_directory', type=str, help='CellRanger Structured FastQ Output Directory')
    args = vars(argparser.parse_args())
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
