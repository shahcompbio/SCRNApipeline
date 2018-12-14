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
import sys
import matplotlib.pyplot as plt
import scanpy.api as sc

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
from interface.singlecellexperiment import SingleCellExperiment

from utils.reporting import Results
from utils.config import Configuration
from utils.export import exportMD, ScaterCode
from utils import plotting

from workflow import PrimaryRun, SecondaryAnalysis

config = Configuration()

def create_workflow():

    workflow = pypeliner.workflow.Workflow()

    bcl_directory = args.get("bcl", None)
    fastq_directories = args.get("fastq", [])
    aggregate = args.get("aggregate_mlibs", list())
    libbase = args.get("lib_base", None)
    combine_assign = args.get("combine",[])
    prefix = args.get("prefix","./")
    output = args.get("out","./")
    recipe = args.get("recipe","basic")

    results = Results(output)
    runner  = PrimaryRun(workflow, prefix, output)

    """
    Aggregating Libraries
    """

    if aggregate != None and len(aggregate) > 0:
        if agg_type == "tenx":
            runner.aggregate_libraries_tenx(aggregate, libbase)
            args["tenx"] = os.path.join(output,"run_{}/outs".format(prefix))
        if agg_type == "scanorama":
            runner.aggregate_libraries_scanorama(aggregate,libbase)
            args["tenx"] = os.path.join(output,"run_{}/outs".format(prefix))

    """
    Setup
    """

    bcls     = runner.set_bcl(bcl_directory)
    fastqs   = runner.set_fastq(fastq_directories)
    workflow = runner.get_workflow()

    tenx_analysis = args.get("tenx", None)
    rdata = args.get("rdata", None)

    secondary_analysis  = SecondaryAnalysis(workflow, prefix, output)
    tenx = TenxAnalysis(tenx_analysis)

    """
    QC
    """

    secondary_analysis.run_scater()
    secondary_analysis.build_sce(tenx)
    secondary_analysis.set_rdata(rdata)

    results.add_analysis(tenx_analysis)
    results.add_workflow(secondary_analysis.rscript)
    results.add_sce(secondary_analysis.sce)

    umi = os.path.join(output,"figures/umi_distribution.png")
    mito = os.path.join(output,"figures/mito_distribution.png")
    ribo = os.path.join(output, "figures/ribo_distribution.png")
    freq = os.path.join(output, "figures/highestExprs.png")
    tech = os.path.join(output, "figures/mean_variance_trend.png")

    results.add_plot(umi,"UMI Distribution")
    results.add_plot(mito,"Mito Distribution")
    results.add_plot(ribo,"Ribo Distribution")
    results.add_plot(freq,"Highest Frequency")
    results.add_plot(tech,"Mean Variance Trend")

    results.add_cellassign_pkl(secondary_analysis.cell_assign_fit)
    results.add_cellassign_raw(secondary_analysis.cell_assign_rdata)

    """
    Cell Assign
    """
    if config.run_cellassign:
        tenx = TenxAnalysis(tenx_analysis)
        if hasattr(config, "rho_matrix"):
            rho_matrix = eval(open(config.rho_matrix,"r").read())
        elif hasattr(config, "organ"):
            sce = SingleCellExperiment.fromRData(secondary_analysis.sce)
            rho_matrix = generate_json(tenx, sce, config.organ)
        else:
            raise AssertionError("Not implemented.")
        secondary_analysis.run_cell_assign(rho_matrix, tenx_analysis, additional=combine_assign)

        path = secondary_analysis.plot_cell_types()
        results.add_plot(path, "Cell Type Frequency")
        path = secondary_analysis.plot_cell_type_by_cluster(tenx_analysis)
        results.add_plot(path, "Cell Type by Cluster")

        path = secondary_analysis.plot_tsne_by_cell_type()
        results.add_plot(path, "TSNE by Cell Type")

        path = secondary_analysis.plot_pca_by_cell_type()
        results.add_plot(path, "TSNE by Cell Type")

    """
    SCVis
    """
    if config.run_scvis:
        secondary_analysis.run_scviz(config.perplexity, 2)
        secondary_analysis.run_scviz(config.perplexity, 10)
        secondary_analysis.run_scviz(config.perplexity, 50)

    """
    Cluster Analysis
    """
    if config.clustering:
        path = secondary_analysis.plot_pca_by_cluster(tenx_analysis, pcs=50)
        results.add_plot(path, "PCA by Cluster")

        path = secondary_analysis.plot_tsne_by_cluster(tenx_analysis, pcs=50)
        results.add_plot(path, "TSNE by Cluster")

        path = secondary_analysis.plot_umap_by_cluster(tenx_analysis, pcs=50)
        results.add_plot(path, "UMAP by Cluster")

        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="PCA", pcs=50)

        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="TSNE", pcs=50)

        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="UMAP", pcs=50)

        template = os.path.join(output,"scvis/5_2/*0.tsv")
        embedding_file = glob.glob(template)[0]
        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="SCVIS", pcs=50, embedding_file=embedding_file)

        template = os.path.join(output,"scvis/5_10/*0.tsv")
        embedding_file = glob.glob(template)[0]
        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="SCVIS", pcs=50, embedding_file=embedding_file)

        template = os.path.join(output,"scvis/5_50/*0.tsv")
        embedding_file = glob.glob(template)[0]
        secondary_analysis.plot_cluster_markers(tenx_analysis, rep="SCVIS", pcs=50, embedding_file=embedding_file)


    if config.plot_scvis:
        template = os.path.join(output,"scvis/5_2/*0.tsv")
        embedding_file = glob.glob(template)[0]
        path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=50)
        path = os.path.join(output, path)
        results.add_plot(path, "SCVis by Cluster (Dim 2)")


        if os.path.exists(secondary_analysis.cell_assign_fit):
            path = secondary_analysis.plot_scvis_by_cell_type(embedding_file, pcs=50)
            results.add_plot(path, "SCVIS (Dim 2) by Cell Type")

        template = os.path.join(output,"scvis/5_10/*0.tsv")
        embedding_file = glob.glob(template)[0]
        path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=50)
        path = os.path.join(output, path)
        results.add_plot(path, "SCVis by Cluster (Dim 10)")


        if os.path.exists(secondary_analysis.cell_assign_fit):
            path = secondary_analysis.plot_scvis_by_cell_type(embedding_file, pcs=50)
            results.add_plot(path, "SCVIS (Dim 10) by Cell Type")

        template = os.path.join(output,"scvis/5_50/*0.tsv")
        embedding_file = glob.glob(template)[0]
        path = secondary_analysis.plot_scvis_by_cluster(tenx_analysis, embedding_file, pcs=50)
        path = os.path.join(output, path)
        results.add_plot(path, "SCVis by Cluster (Dim 50)")


        if os.path.exists(secondary_analysis.cell_assign_fit):
            path = secondary_analysis.plot_scvis_by_cell_type(embedding_file, pcs=50)
            results.add_plot(path, "SCVIS (Dim 50) by Cell Type")

    """
    Gene Level
    """
    ## Table of Data

    ## Differential  Genes by Cluster

    ## Differential Genes by Cell Type

    """
    Clone Align
    """

    """
    Reporting
    """
    if config.report:
        workflow.transform (
            name = "{}_markdown".format(prefix),
            func = exportMD,
            args = (
                results,
            )
        )



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
