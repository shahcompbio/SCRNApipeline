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
import yaml

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
from interface.genemarkermatrix import GeneMarkerMatrix, generate_json

from utils.reporting import Results
from utils.config import *
from utils.export import exportMD, ScaterCode
from utils import plotting

from workflow import PrimaryRun, SecondaryAnalysis

import sys


def yaml_configuration():
    yaml_file = args.get("yaml",None)
    if yaml_file is not None:
        with open(yaml_file, "r") as f:
            doc = yaml.load(f)
            for var in doc:
                print(var)


def create_workflow():

    workflow = pypeliner.workflow.Workflow()
    yaml_configuration()

    bcl_directory = args.get("bcl", None)
    fastq_directories = args.get("fastq", [])
    aggregate = args.get("aggregate_mlibs", list())
    libbase = args.get("lib_base", None)
    combine_assign = args.get("combine",[])
    prefix = args.get("prefix","./")
    output = args.get("out","./")
    recipe = args.get("recipe","basic")

    results = Results(output, prefix)

    runner   = PrimaryRun(workflow, prefix, output)

    if aggregate != None and len(aggregate) > 0:
        runner.aggregate_libraries(aggregate, libbase)
        workflow = runner.get_workflow()
        return workflow

    bcls     = runner.set_bcl(bcl_directory)
    fastqs   = runner.set_fastq(fastq_directories)
    workflow = runner.get_workflow()

    tenx_analysis = args.get("tenx", None)
    rdata = args.get("rdata", None)

    secondary_analysis  = SecondaryAnalysis(workflow, prefix, output)
    tenx = TenxAnalysis(tenx_analysis)

    secondary_analysis.run_scater()
    secondary_analysis.build_sce(tenx)
    secondary_analysis.set_rdata(rdata)

    results.add_workflow(secondary_analysis.rscript)
    results.add_filtered_sce(secondary_analysis.sce_filtered)
    results.add_final_sce(secondary_analysis.sce_final)


    umi = os.path.join(output,"umi_distribution.png")
    mito = os.path.join(output,"mito_distribution.png")
    ribo = os.path.join(output, "ribo_distribution.png")
    freq = os.path.join(output, "highestExprs.png")

    results.add_plot(umi,"UMI Distribution")
    results.add_plot(mito,"Mito Distribution")
    results.add_plot(ribo,"Ribo Distribution")
    results.add_plot(freq,"Highest Frequency")

    # rho_matrix = generate_json(this_data)

    # analysis.run_cell_assign(rho_matrix, tenx_analysis, additional=combine_assign)
    #
    # results.add_cellassign_pkl(analysis.cell_assign_fit)
    # results.add_cellassign_raw(analysis.cell_assign_rdata)

    # perplexity = 5
    # component = 2
    # analysis.run_scviz(perplexity, 2)
    # analysis.run_scviz(perplexity, 10)
    # analysis.run_scviz(perplexity, 50)
    # tenx = TenxAnalysis(tenx_analysis)

    path = secondary_analysis.plot_tsne_by_cluster(tenx_analysis)
    path = secondary_analysis.plot_tsne_by_cluster(tenx_analysis, raw=True)
    template = os.path.join(output,"5_2_raw/*0.tsv")
    embedding_file = glob.glob(template)[0]
    path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=2, raw=True)

    template = os.path.join(output,"5_10_raw/*0.tsv")
    embedding_file = glob.glob(template)[0]
    path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=10, raw=True)

    template = os.path.join(output,"5_50_raw/*0.tsv")
    embedding_file = glob.glob(template)[0]
    path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=50, raw=True)


    # path = secondary_analysis.plot_tsne_by_cluster(tenx_analysis, rep="UMAP", raw=True)
    # path = secondary_analysis.plot_tsne_by_cluster(tenx_analysis, raw=True)

    # results.add_plot(path, "TSNE by Cluster")
    # path = analysis.plot_tsne_by_cell_type()
    # results.add_plot(path, "TSNE by Cell Type")
    # path = analysis.plot_cell_types()
    # results.add_plot(path, "Cell Type Frequency")
    # path = analysis.plot_cell_type_by_cluster()
    # results.add_plot(path, "Cell Type by Cluster")
    # results.barcode_to_celltype()
    # workflow.transform (
    #     name = "{}_markdown".format(prefix),
    #     func = exportMD,
    #     args = (
    #         results,
    #     )
    # )

    workflow = secondary_analysis.get_workflow()
    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('--bcls', type=str, help='BaseCalls Illumina Directory')
    argparser.add_argument('--fastqs', type=str, nargs="+", help='CellRanger Structured FastQ Output Directory')
    argparser.add_argument('--tenx', type=str, help='Output Directory From Cell Ranger mkfastq - or any folder with *_bc_gene_matrices')
    argparser.add_argument('--rdata', type=str, help='Serialized Single Cell Experiment From R')
    argparser.add_argument('--out', type=str, help="Base directory for output")
    argparser.add_argument("--prefix", type=str, help="Analysis prefix")
    argparser.add_argument("--aggregate-mlibs", nargs='+', type=str, help="Library prefixes to aggregate.")
    argparser.add_argument("--lib-base", type=str, help="Directory containing libraries to aggregate by prefix (must provide aggregate-mlibs)")
    argparser.add_argument("--combine", nargs='+', type=str, help="Library prefixes to aggregate.")
    argparser.add_argument("--yaml", type=str, help="Configuration settings for pipeline.")
    argparser.add_argument("--recipe", type=str, help="Normalization and filtering recipe.")

    parsed_args = argparser.parse_args()

    args = vars(parsed_args)
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
