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

import argparse
import glob
import os

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

from software.cellranger import CellRanger
from software.tenx import TenX
from software.cellassign import CellAssign
from software.clonealign import CloneAlign
from software.scviz import SCViz
from software.fastqc import FastQC

from interface.binarybasecall import BinaryBaseCall
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis
from interface.qcreport import QCReport
from interface.genemarkermatrix import GeneMarkerMatrix

from utils.reporting import HTMLResults
from utils.config import *
from utils.export import exportRMD



def create_run_workflow(subprefix, output, bcl_directory, fastq_directory):

    workflow = pypeliner.workflow.Workflow()

    bcl = None
    fastq = None
    tenx = None
    qcreport = QCReport(output)
    qcscript = qcreport.generate_script()

    # try:
    #     os.makedirs(fastqc_output)
    # except Exception as e:
    #     pass
    #
    # if bcl_directory:
    #     bcl = BinaryBaseCall(bcl_directory)
    #     workflow.transform (
    #         name = "{}_cellranger_mkfastq".format(subprefix),
    #         func = CellRanger.mkfastq,
    #         ret = pypeliner.managed.TempOutputObj("fastq_object_{}".format(subprefix)),
    #         args = (
    #             bcl,
    #         )
    #     )
    #

    if fastq_directory != None:
        fastq = FastQDirectory(fastq_directory, subprefix, output)
    else:
        fastq = pypeliner.managed.TempInputObj("fastq_object_{}".format(subprefix))
    if not fastq.has_qc():
        workflow.transform (
            name = "{}_fastqc".format(subprefix),
            func = FastQC.run,
            args = (
                fastq,
            )
        )

    if not fastq.check_status():
        workflow.transform (
            name = "{}_cellranger_counts".format(subprefix),
            func = CellRanger.count,
            ret = pypeliner.managed.TempOutputObj("tenx_analysis_{}".format(subprefix)),
            args = (
                fastq,
            )
        )
    if not tenx:
        tenx = TenxAnalysis(fastq.results)
    else:
        tenx = pypeliner.managed.TempInputObj("tenx_analysis_{}".format(subprefix))
    print("Running...")
    # workflow.transform (
    #     name = "{}_read10xcountsFiltered".format(subprefix),
    #     func = TenX.read10xCountsFiltered,
    #     args = (
    #         tenx,
    #         pypeliner.managed.OutputFile("sce_filtered.rdata")
    #     )
    # )
    # #
    # workflow.commandline (
    #     name = "{}_scater".format(subprefix),
    #     args = ("Rscript", qcscript, pypeliner.managed.InputFile("sce_filtered.rdata"), pypeliner.managed.OutputFile("sce_final.rdata")),
    # )

    if rho_matrix is not None and not os.path.exists("cell_assign_fit.pkl"):
        workflow.transform (
            name = "{}_cellassign".format(subprefix),
            func = CellAssign.run_em,
            args = (
                "sce_final.rdata",
                pypeliner.managed.OutputFile("cell_assign_fit.pkl"),
                subprefix
            )
        )

    # if copy_number_data is not None:
    #     workflow.transform (
    #         name = "{}_clonealign".format(subprefix),
    #         func = CloneAlign.run,
    #         args = (
    #             pypeliner.managed.InputFile("sce_final.rdata"),
    #             copy_number_data, #TODO  move copy number loading code
    #             pypeliner.managed.OutputFile("clone_align_fit.rdata")
    #         )
    #     )
    # if scviz_embedding is not None:
    #     workflow.transform (
    #         name = "{}_scviz".format(subprefix),
    #         func = SCViz.run,
    #         args = (
    #             pypeliner.managed.InputFile("sce_final.rdata"),
    #             scviz_embedding,
    #             output
    #         )
    #     )
    # else:
    # print("Running SC")
    # workflow.transform (
    #     name = "{}_scviz".format(subprefix),
    #     func = SCViz.train,
    #     args = (
    #         pypeliner.managed.InputFile("sce_final.rdata"),
    #         tenx,
    #         pypeliner.managed.InputFile("cell_assign_fit.pkl")
    #     )
    # )

    workflow.transform (
        name = "{}_save_rmd".format(subprefix),
        func = exportRMD,
        args = (
            output,
            subprefix,
            fastq,
            pypeliner.managed.InputFile("sce_final.rdata")
        )
    )

    return workflow


def create_workflow():

    workflow = pypeliner.workflow.Workflow()

    bcl_directory = args.get("bcl", None)
    fastq_directory = args.get("fastq", None)
    tenx_analysis = args.get("tenx", None)
    rdata = args.get("rdata", None)
    output = args.get("out","./")
    prefix = args.get("prefix","./")

    try:
        os.makedirs(os.path.join(output,"fastqc"))
        print("FastQC directory created.")
    except OSError:
        print("FastQC directory already created.")

    if bcl_directory:
        bcl_dirs = glob.glob(os.path.join(bcl_directory, "**/*.csv"), recursive=True)
        for i, bcl_dir in enumerate(bcl_dirs):
            subprefix = "{}_{}".format(prefix,i)
            workflow.subworkflow(
                name = subprefix,
                func = create_run_workflow,
                args = (
                    subprefix,
                    output,
                    os.path.split(bcl_directory)[0],
                    None
                )
            )
    else:
        print(fastq_directory)
        fastq_dirs = glob.glob(os.path.join(fastq_directory, "**/*.csv"),recursive=True)
        for i, fastq_dir in enumerate(fastq_dirs):
            subprefix = "{}_{}".format(prefix,i)
            workflow.subworkflow(
                name = subprefix,
                func = create_run_workflow,
                args = (
                    subprefix,
                    output,
                    None,
                    os.path.split(fastq_dir)[0]
                )
            )

    return workflow

if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--bcl', type=str, help='BaseCalls Illumina Directory')
    argparser.add_argument('--fastq', type=str, help='CellRanger Structured FastQ Output Directory')
    argparser.add_argument('--tenx', type=str, help='Output Directory From Cell Ranger mkfastq - or any folder with *_bc_gene_matrices')
    argparser.add_argument('--rdata', type=str, help='Serialized Single Cell Experiment From R')
    argparser.add_argument('--out', type=str, help="Base directory for output")
    argparser.add_argument("--prefix", type=str, help="Analysis prefix")

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
