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
from interface.genemarkermatrix import GeneMarkerMatrix

from utils.reporting import HTMLResults
from utils.config import *
from utils.export import exportRMD, ScaterCode, cat
from utils import plotting



def create_run_workflow(prefix, output, bcl, fastqs, aggregate_csv):

    workflow = pypeliner.workflow.Workflow()

    # if bcl is not None:
    #     workflow.transform (
    #         name = "{}_cellranger_mkfastq".format(prefix),
    #         func = CellRanger.mkfastq,
    #         args = (
    #             bcl,
    #         )
    #     )
    #     fastq = FastQDirectory(bcl.out(), subprefix, output)

    # for fastq in fastqs:
    #     if not fastq.has_qc():
    #         workflow.transform (
    #             name = "{}_fastqc".format(prefix),
    #             func = FastQC.run,
    #             args = (
    #                 fastq,
    #             )
    #         )
    
    if len(fastqs) > 0:
        workflow.transform (
            name = "{}_mro_build".format(prefix),
            
        )
    workflow.transform (
        name = "{}_cellranger_counts".format(prefix),
        func = CellRanger.count,
        args = (
            fastqs,
        )
    )

    # if aggregate_csv is not None:
    #     print("Running Aggr...")
    #     workflow.transform (
    #         name = "{}_aggr".format(prefix),
    #         func = CellRanger.aggr,
    #         args = (
    #             aggregate_csv,
    #             prefix
    #         )
    #     )
    #     cellranger_output = os.path.join(output, prefix, "outs")
    #     analysis = TenxAnalysis(cellranger_output)
    # print("Fastq",fastqs[0].results)
    # fastq = FastQDirectory(fastq_directory, prefix, output)

    # workflow.subworkflow (
    #     name = "{}_analysis".format(prefix),
    #     func = create_analysis_workflow,
    #     args = (
    #         prefix,
    #         analysis,
    #         fastqs,
    #         output
    #     )
    # )

    return workflow


def create_analysis_workflow(prefix, analysis, fastqs, output):

    workflow = pypeliner.workflow.Workflow()

    scater_workflow = ScaterCode(output)
    rscript = scater_workflow.generate_script()

    sce_filtered = os.path.join(output,"sce_filtered.rdata")
    sce_final = os.path.join(output,"sce_final.rdata")
    cell_assign_fit = os.path.join(output, "cell_assign_fit.pkl")
    clone_align_fit = os.path.join(output, "clone_align_fit.rdata")

    # if not os.path.exists(sce_filtered):
    #     workflow.transform (
    #         name = "{}_read10xcounts".format(prefix),
    #         func = TenX.read10xCountsFiltered,
    #         args = (
    #             analysis,
    #             pypeliner.managed.OutputFile(sce_filtered)
    #         )
    #     )

    # print("Running RSCript")
    # workflow.commandline (
    #     name = "{}_scater_workflow".format(prefix),
    #     args = ("Rscript", rscript, pypeliner.managed.InputFile(sce_filtered), pypeliner.managed.OutputFile(sce_final)),
    # )

    # if rho_matrix is not None:
    #     workflow.transform (
    #         name = "{}_cellassign".format(prefix),
    #         func = CellAssign.run_em,
    #         args = (
    #             pypeliner.managed.InputFile(sce_final),
    #             cell_assign_fit,
    #             prefix
    #         )
    #     )

    # workflow.transform (
    #     name = "celltypes",
    #     func = plotting.celltypes,
    #     args = (
    #         pypeliner.managed.InputFile(sce_final),
    #         pypeliner.managed.InputFile(cell_assign_fit),
    #     )
    # )

    workflow.transform (
        name = "tsne_by_cluster",
        func = plotting.tsne_by_cluster,
        args = (
            pypeliner.managed.InputFile(sce_final),
            pypeliner.managed.InputFile(cell_assign_fit),
        )
    )

    workflow.transform (
        name = "tsne_by_cell_type",
        func = plotting.tsne_by_cell_type,
        args = (
            pypeliner.managed.InputFile(sce_final),
            pypeliner.managed.InputFile(cell_assign_fit),
        )
    )

    workflow.transform (
        name = "cell_type_by_cluster",
        func = plotting.cell_type_by_cluster,
        args = (
            pypeliner.managed.InputFile(sce_final),
            pypeliner.managed.InputFile(cell_assign_fit),
        )
    )

    # workflow.transform (
    #     name = "{}_umap".format(prefix),
    #     func = plotting.umap,
    #     args = (
    #         pypeliner.managed.InputFile(sce_final),
    #         analysis,
    #         pypeliner.managed.InputFile(cell_assign_fit),
    #         "umap.png"
    #     )
    # )


    # if copy_number_data is not None:
    #     workflow.transform (
    #         name = "{}_clonealign".format(prefix),
    #         func = CloneAlign.run,
    #         args = (
    #             pypeliner.managed.InputFile(sce_final),
    #             copy_number_data, #TODO  move copy number loading code
    #             pypeliner.managed.OutputFile(clone_align_fit)
    #         )
    #     )
    #
    #
    # if scviz_embedding is not None:
    #     workflow.transform (
    #         name = "{}_scviz".format(prefix),
    #         func = SCViz.run,
    #         args = (
    #             pypeliner.managed.InputFile(sce_final),
    #             scviz_embedding,
    #             output
    #         )
    #     )
    # else:
    #     workflow.transform (
    #         name = "{}_scviz".format(prefix),
    #         func = SCViz.train,
    #         args = (
    #             pypeliner.managed.InputFile(sce_final),
    #             analysis,
    #             pypeliner.managed.InputFile(cell_assign_fit)
    #         )
    #     )

    output = args.get("out","./")
    workflow.transform (
        name = "{}_save_rmd".format(prefix),
        func = exportRMD,
        args = (
            fastqs[0],
            analysis,
            scater_workflow,
            prefix,
            sce_final,
            output
        )
    )

    return workflow


def create_workflow():

    workflow = pypeliner.workflow.Workflow()

    bcl_directory = args.get("bcl", None)
    fastq_directories = args.get("fastq", [])
    tenx_analysis = args.get("tenx", None)
    rdata = args.get("rdata", None)
    output = args.get("out","./")
    prefix = args.get("prefix","./")
    aggregate = args.get("aggregate-mlibs", list())

    try:
        os.makedirs(os.path.join(output,"fastqc"))
        print("FastQC directory created.")
    except OSError:
        print("FastQC directory already created.")

    if bcl_directory:
        bcl_subdirectories = glob.glob(os.path.join(bcl_directory, "**/*.csv"), recursive=True)
        for i, bcl_subdirectory in enumerate(bcl_subdirectories):
            subprefix = "{}_{}".format(prefix,i)
            bcl = BinaryBaseCall(os.path.split(bcl_subdirectory)[0])
            workflow.subworkflow(
                name = subprefix,
                func = create_run_workflow,
                args = (
                    subprefix,
                    output,
                    bcl,
                    None,
                    None
                )
            )

    elif len(fastq_directories) > 0:
        fastqs = []
        for fastq_directory in fastq_directories:
            fastqs.append(FastQDirectory(fastq_directory, prefix, output))
        workflow.subworkflow(
            name = prefix,
            func = create_run_workflow,
            args = (
                prefix,
                output,
                None,
                fastqs,
                None
            )
        )

    elif combine:
        raise AssertionError("Not Implemented.")

    elif aggregate:
        aggregated_output = os.path.join(output,prefix)
        libraries = os.path.join(aggregated_output,"libraries.csv")
        csv = open(libraries,"w")
        csv.write("library_id,molecule_h5\n")
        for library in aggregate:
            lib_output = os.path.join(output,library)
            partial_path = os.path.join(lib_output,"*/outs/molecule_info.h5")
            for molecule_h5 in glob.glob(partial_path):
                if os.path.exists(molecule_h5):
                    csv.write("{prefix},{path}\n".format(prefix=library,path=molecule_h5))
                else:
                    print("Not Found", molecule_h5)
        csv.close()
        workflow.subworkflow(
            name = "aggr_{}".format(prefix),
            func = create_run_workflow,
            args = (
                prefix,
                aggregated_output,
                None,
                None,
                libraries
            )
        )

    return workflow



if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--bcl', type=str, help='BaseCalls Illumina Directory')
    argparser.add_argument('--fastq', type=str, nargs="+", help='CellRanger Structured FastQ Output Directory')
    argparser.add_argument('--tenx', type=str, help='Output Directory From Cell Ranger mkfastq - or any folder with *_bc_gene_matrices')
    argparser.add_argument('--rdata', type=str, help='Serialized Single Cell Experiment From R')
    argparser.add_argument('--out', type=str, help="Base directory for output")
    argparser.add_argument("--prefix", type=str, help="Analysis prefix")
    argparser.add_argument("--aggregate-mlibs", nargs='+', type=str, help="Library prefixes to aggregate.")
    # argparser.add_argument("--aggregate-slibs", nargs='+', type=str, help="Sample prefixes to aggregate.")

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
