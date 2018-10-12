"""

Single Cell RNA-Seq Pipeline

Pipeline for running single cell rna-seq experiments.
This is primarily built on Cell Ranger with additionaly analysis
from CellAssign, CloneAlign, and SCViz tools.
The workflow is inferred based on the inclusion (or omission) of command line arguments.

Example:
    Running pipeline starting from binary base call directory (BCL) to report generation.
    The top level BCL directory must include a single csv worksheet with minimum columns
    for Lane, Sample, Index (10x index used for library construction).

        $ python3 pipeline.py --bcl tests/cellranger-tiny-bcl-1.2.0/

    The FastsQ directory can be the output of `cellranger mkfastq` or a directory
    where fastq files are named '{sample}_S{sample_num}_L00{lane}_{R{read} || I1}_001'.
    If this argument is omitted, but the `bcl` argument is included,
    the fastq path will be inferred.

        $ python3 pipeline.py --fastq tests/tiny-fastqs-mk/outs/fastq_path/

    The tenx analysis folder can be the output of `cellranger count`
    or a directory that includes a {filtered || raw}_gene_bc_matrices folder.
    If this argument is omitted, but bcl or fastq is included, the path will be inferred.
    If this directory includes Cell Ranger analysis,
    this will be included in the report generation.

        $ python3 pipeline.py --tenx tests/tiny-fastqs-count/outs/

    Single Cell Experiment objects can be created from tenx analysis folders
    or loaded from serialized RData objects.
    You can load these and run the pipeline downstream analysis
    starting from these serialized objects.

        $ python3 pipeline.py --rdata tests/example_sce.RData


"""



import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from software.cellranger import CellRanger
from software.tenx import TenX
from software.clonealign import CloneAlign
from software.cellassign import CellAssign
from software.scviz import SCViz

from interface.binarybasecall import BinaryBaseCall
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis

from utils.reporting import HTMLResults

import argparse

def create_workflow():
    """
    Generates tasks as Pypeliner workflow based on input arguments.
    The workflow will start from the most raw input provided and override
    any downstream tasks with subsequently provided input arguments.
    Parellelization is performed over provided samplesheets and replication
    within samples.

    Args:
        None

    Yields:
        Pypeliner workflow object.
    """

    bcl_directory = args.get("bcl", None)
    fastq_directory = args.get("fastq", None)
    tenx_analysis = args.get("tenx", None)
    rdata = args.get("rdata", None)

    bcl = BinaryBaseCall(bcl_directory)

    workflow = pypeliner.workflow.Workflow()

    if bcl_directory:
        workflow.transform (
            name = "bcl_to_fastq",
            func = CellRanger.mkfastq,
            ret = pypeliner.managed.TempOutputObj("fastq_object"),
            args = (
                bcl_object,
            )
        )

    if bcl_directory != None or fastq_directory != None:
        if fastq_directory != None:
            fastq = FastQDirectory(fastq_directory)
        else:
            fastq = pypeliner.managed.TempInputObj("fastq_object")
        workflow.transform (
            name = "fastqc",
            func = FastQC.run,
            ret = pypeliner.managed.TempOutputObj("qc_report")
            args = (
                fastq,
            )
        )

    if bcl_directory != None or fastq_directory != None:
        if fastq_directory != None:
            fastq = FastQDirectory(fastq_directory)
        else:
            fastq = pypeliner.managed.TempInputObj("fastq_object")
        workflow.transform (
            name = "fastq_counts",
            func = CellRanger.count,
            ret = pypeliner.managed.TempOutputObj("tenx_analysis"),
            args = (
                fastq,
            )
        )

    tenx = None
    if tenx_analysis != None and rdata == None:
        tenx = TenxAnalysis(tenx_analysis)
    elif tenx_analysis == None and rdata == None:
        tenx = pypeliner.managed.TempInputObj("tenx_analysis")
    if tenx != None:
        workflow.transform (
            name = "tenx_read10xcounts",
            func = TenX.read10xCounts,
            ret = pypeliner.managed.TempOutputObj("single_cell_experiment"),
            args = (
                tenx,
            )
        )

    if rdata != None:
        single_cell_experiment = TenxAnalysis.from_rdata(rdata)
    else:
        single_cell_experiment = pypeliner.managed.TempInputObj("single_cell_experiment")



    #
    # workflow.transform (
    #     name = "tenx_barcoderanks",
    #     func = TenX.barcodeRanks,
    #     args = (
    #         pypeliner.managed.TempInputObj("tenx_analysis"),
    #     )
    # )
    #
    # workflow.transform (
    #     name = "tenx_emptydrops",
    #     func = TenX.emptyDrops,
    #     args = (
    #         pypeliner.managed.TempInputObj("tenx_analysis"),
    #     )
    # )

    workflow.transform (
        name = "clonealign",
        func = CloneAlign.run,
        ret = pypeliner.managed.TempOutputObj("clone_align_fit"),
        args = (
            single_cell_experiment,
        )
    )

    # """
    # workflow.transform (
    #     name = "cellasign",
    #     func = CellAssign.run_em,
    #     ret = pypeliner.managed.TempOutputObj("cell_assignments"),
    #     args = (
    #         single_cell_experiment,
    #     )
    # )
    #
    # workflow.transform (
    #     name = "scviz",
    #     func = SCViz.run,
    #     ret = pypeliner.managed.TempOutputObj("scviz_dim_reduction"),
    #     args = (
    #         single_cell_experiment,
    #     )
    # )
    #
    # workflow.transform (
    #     name = "html_output",
    #     func = HTMLResults.generate,
    #     args = (
    #         pypeliner.managed.TempInputObj("fastq_object")
    #         pypeliner.managed.TempInputObj("ten_analysis"),
    #         pypeliner.managed.TempInputObj("single_cell_experiment"),
    #         pypeliner.managed.TempInputObj("clone_align_fit"),
    #         pypeliner.managed.TempInputObj("cell_assignments"),
    #         pypeliner.managed.TempInputObj("scviz_dim_reduction"),
    #     )
    # )
    """

    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--bcl', type=str, help='BaseCalls Illumina Directory')
    argparser.add_argument('--fastq', type=str, help='CellRanger Structured FastQ Output Directory')
    argparser.add_argument('--tenx', type=str, help='Output Directory From Cell Ranger mkfastq - or any folder with *_bc_gene_matrices')
    argparser.add_argument('--rdata', type=str, help='Serialized Single Cell Experiment From R')
    argparser.add_argument('--nosecondary', type=str, help="Disable all downstream analysis")

    args = vars(argparser.parse_args())
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
