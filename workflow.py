import os
import glob

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
from interface.singlecellexperiment import SingleCellExperiment
from utils.export import exportRMD, ScaterCode

from utils import plotting, combine

import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

class PrimaryRun(object):

    def __init__(self, workflow, prefix, output):
        self.prefix = prefix
        self.subprefixes = []
        self.workflow = workflow
        self.output = output
        try:
            os.makedirs(os.path.join(output,"fastqc"))
        except OSError:
            pass

    def generate_subprefix(self, directory):
        return "generic"

    def set_bcl(self, bcl_directory):
        if bcl_directory is not None:
            self.bcl = bcl_directory
            self.bcl_directories = []
            subdirectories = glob.glob(os.path.join(bcl_directory, "**/*.csv"), recursive=True)
            fastq_directories = []
            for subdirectory in subdirectories:
                subprefix = generate_subprefix(directory)
                bcl = BinaryBaseCall(os.path.split(subdirectory)[0])
                self.workflow.transform (
                    name = "{}_cellranger_mkfastq".format(subprefix),
                    func = CellRanger.mkfastq,
                    args = (
                        bcl,
                    )
                )
                fastq_directories.append(bcl.out())
            self.set_fastq(fastq_directories)
        else:
            print("No BCL to run.")

    def set_fastq(self, fastq_directories):
        self.fastq_directories = fastq_directories
        self.fastqs = []
        for fastq_directory in self.fastq_directories:
            fastq = FastQDirectory(fastq_directory, prefix, output)
            #if not fastq.has_qc():
            if False:
                self.workflow.transform (
                    name = "{}_fastqc".format(self.prefix),
                    func = FastQC.run,
                    args = (
                        fastq,
                    )
                )
            if fastq.check_status():
                self.fastqs.append(fastq)
        if len(self.fastqs) > 0:
            self.workflow.transform (
                name = "{}_cellranger_counts".format(self.prefix),
                func = CellRanger.count,
                args = (
                    self.fastqs,
                )
            )
        else:
            print("No FastQ to run.")

    def _generate_csv(self, libbase):
        csv = open(self.libraries_csv,"w")
        csv.write("library_id,molecule_h5\n")
        for library in self.libraries:
            print(library)
            lib_output = os.path.join(libbase,library)
            print(lib_output)
            molecule_h5 = os.path.join(lib_output,"outs/molecule_info.h5")
            if os.path.exists(molecule_h5):
                print("Writing")
                csv.write("{prefix},{path}\n".format(prefix=library,path=molecule_h5))
            else:
                print("Not Found", molecule_h5)
        csv.close()

    def aggregate_libraries_tenx(self, libraries, libbase):
        if len(libraries) > 0:
            self.libraries = libraries
            self.libraries_csv = os.path.join(self.output, "libraries.csv")
            self._generate_csv(libbase)
            self.workflow.transform (
                name = "{}_aggr".format(self.prefix),
                func = CellRanger.aggr,
                args = (
                    self.libraries_csv,
                    self.prefix
                )
            )

    def get_workflow(self):
        return self.workflow

class SecondaryAnalysis(object):

    def __init__(self, workflow, prefix, output):
        self.workflow = workflow
        self.prefix = prefix
        self.output = output
        self.raw_matrices_sce = os.path.join(self.output,"rdata/raw_matrices_sce.rdata")
        self.filtered_matrices_sce = os.path.join(self.output,"rdata/filtered_matrices_sce.rdata")
        self.sce = os.path.join(self.output,"rdata/sce.rdata")
        self._pca = os.path.join(self.output,"rdata/pca_{}.rdata")
        self._raw_pca = os.path.join(self.output,"rdata/raw_pca_{}.rdata")
        self.cell_assign_fit = os.path.join(self.output, "rdata/cell_assign_fit.pkl")
        self.cell_assign_rdata = os.path.join(self.output, "rdata/cell_assign_fit.rdata")
        self.clone_align_fit = os.path.join(self.output, "rdata/clone_align_fit.rdata")
        self.cell_types_fig = os.path.join(self.output,"figures/cell_types.png")


    def set_rdata(self, rdata):
        if rdata is not None:
            self.rdata = rdata

    def build_sce(self,tenx):
        self.workflow.transform (
            name = "read10xcounts",
            func = TenX.read10xCountsFiltered,
            args = (
                tenx,
                pypeliner.managed.OutputFile(self.filtered_matrices_sce)
            )
        )
        self.rdata = self.filtered_matrices_sce

    def save_raw_sce(self,tenx):
        self.workflow.transform (
            name = "read10xcounts_raw",
            func = TenX.read10xCountsRaw,
            args = (
                tenx,
                pypeliner.managed.OutputFile(self.raw_matrices_sce)
            )
        )

    def run_scater(self):
        self.scater_workflow = ScaterCode(self.output)
        self.rscript = self.scater_workflow.generate_script()
        if not os.path.exists(self.sce):
            self.workflow.commandline (
                name = "{}_scater_workflow".format(self.prefix),
                args = ("Rscript", self.rscript, pypeliner.managed.InputFile(self.filtered_matrices_sce), pypeliner.managed.OutputFile(self.sce)),
            )
        self.rdata = self.sce

    def run_cell_assign(self, rho_matrix, tenx, additional=None):
        self.workflow.transform (
            name = "{}_cellassign".format(self.prefix),
            func = CellAssign.run_em,
            args = (
                tenx,
                pypeliner.managed.InputFile(self.sce),
                pypeliner.managed.OutputFile(self.cell_assign_fit),
                self.prefix,
                rho_matrix,
                additional
            )
        )

    def run_clone_align(self, copy_number_data):
        workflow.transform (
            name = "{}_clonealign".format(self.prefix),
            func = CloneAlign.run,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                copy_number_data,
                pypeliner.managed.OutputFile(clone_align_fit)
            )
        )


    def plot_cell_types(self):
        if not os.path.exists("figures/cell_types.png"):
            self.workflow.transform (
                name = "celltypes",
                func = plotting.celltypes,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    pypeliner.managed.InputFile(self.cell_assign_fit),
                    self.prefix,
                    self.output
                )
            )
        return "figures/cell_types.png"

    def plot_tsne_by_cluster(self, tenx, pcs=2):
        if not os.path.exists("figures/tsne_by_cluster.png"):
            self.workflow.transform (
                name = "tsne_by_cluster",
                func = plotting.tsne_by_cluster,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/tsne_by_cluster.png"

    def tsne_by_cluster_markers(self, tenx, pcs=50):
        if not os.path.exists("figures/rank_genes_groups_leiden_TSNE_2.png"):
            self.workflow.transform (
                name = "tsne_by_cluster_markers",
                func = plotting.tsne_by_cluster_markers,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/rank_genes_groups_leiden_TSNE_2.png"

    def pca_by_cluster_markers(self, tenx, pcs=50):
        if not os.path.exists("figures/rank_genes_groups_leiden_PCA_{}.png".format(pcs)):
            self.workflow.transform (
                name = "pca_by_cluster_markers",
                func = plotting.pca_by_cluster_markers,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/rank_genes_groups_leiden_PCA_{}.png".format(pcs)

    def plot_cluster_markers(self, tenx, rep="TSNE", pcs=50, embedding_file = None):
        if embedding_file is not None:
            if "5_2" in embedding_file:
                suffix = "5_2"
            if "5_10" in embedding_file:
                suffix = "5_10"
            if "5_50" in embedding_file:
                suffix= "5_50"
            name = "cluster_markers_{}_{}".format(rep,suffix)
        else:
            name = "cluster_markers_{}".format(rep)
        self.workflow.transform (
            name = name,
            func = plotting.cluster_markers,
            args = (
                pypeliner.managed.InputFile(self.sce),
                tenx,
                rep,
                pcs,
                embedding_file,
                name.replace("cluster_markers_","")
            )
        )
        return "figures/expression"


    def umap_by_cluster_markers(self, tenx, pcs=50):
        if not os.path.exists("figures/rank_genes_groups_leiden_UMAP_2.png"):
            self.workflow.transform (
                name = "umap_by_cluster_markers",
                func = plotting.umap_by_cluster_markers,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/rank_genes_groups_leiden_UMAP_2.png"

    def scvis_by_cluster_markers(self, tenx, embedding_file, pcs=2):
        if not os.path.exists("figures/rank_genes_groups_leiden_SCVIS_{}.png".format(pcs)):
            self.workflow.transform (
                name = "scvis_by_cluster_markers_{}".format(pcs),
                func = plotting.scvis_by_cluster_markers,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs,
                    embedding_file
                )
            )
        return "figures/rank_genes_groups_leiden_SCVIS_{}.png".format(pcs)

    def plot_pca_by_cluster(self, tenx, pcs=2):
        if not os.path.exists("figures/pca_by_cluster.png"):
            self.workflow.transform (
                name = "pca_by_cluster",
                func = plotting.pca_by_cluster,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/pca_by_cluster.png"

    def plot_scvis_by_cluster(self, tenx, embedding_file, pcs=2):
        if embedding_file is not None:
            if "5_2" in embedding_file:
                suffix = "5_2"
            if "5_10" in embedding_file:
                suffix = "5_10"
            if "5_50" in embedding_file:
                suffix= "5_50"
            prefix = self.prefix + "_{}".format(suffix)
        else:
            prefix = self.prefix
        if not os.path.exists("figures/scvis_by_cluster_{}.png".format(prefix)):
            self.workflow.transform (
                name = "scvis_by_cluster_{}".format(suffix),
                func = plotting.scvis_by_cluster,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    prefix,
                    pcs,
                    embedding_file,
                )
            )
        return "figures/scvis_by_cluster_{}.png".format(prefix)

    def plot_tsne_by_cell_type(self):
        if not os.path.exists("figures/tsne_by_cell_type.png"):
            self.workflow.transform (
                name = "tsne_by_cell_type",
                func = plotting.tsne_by_cell_type,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    pypeliner.managed.InputFile(self.cell_assign_fit),
                    self.prefix,
                )
            )
        return "figures/tsne_by_cell_type.png"

    def plot_pca_by_cell_type(self):
        if not os.path.exists("figures/pca_by_cell_type.png"):
            self.workflow.transform (
                name = "pca_by_cell_type",
                func = plotting.pca_by_cell_type,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    pypeliner.managed.InputFile(self.cell_assign_fit),
                    self.prefix,
                )
            )
        return "figures/pca_by_cell_type.png"


    def plot_scvis_by_cell_type(self, embedding_file, pcs=2):
        if "5_2" in embedding_file:
            suffix = "5_2"
        if "5_10" in embedding_file:
            suffix = "5_10"
        if "5_50" in embedding_file:
            suffix= "5_50"
        prefix = self.prefix + "_{}".format(suffix)
        if not os.path.exists("figures/scvis_by_cell_type_{}.png".format(prefix)):
            self.workflow.transform (
                name = "scvis_by_cell_type_{}".format(suffix),
                func = plotting.scvis_by_cell_type,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    pypeliner.managed.InputFile(self.cell_assign_fit),
                    prefix,
                    embedding_file
                )
            )
        return "figures/scvis_by_cell_type_{}.png".format(prefix)

    def plot_cell_type_by_cluster(self, tenx):
        if not os.path.exists("figures/cell_type_by_cluster.png"):
            self.workflow.transform (
                name = "cell_type_by_cluster",
                func = plotting.cell_type_by_cluster,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    pypeliner.managed.InputFile(self.cell_assign_fit),
                    tenx,
                    self.prefix
                )
            )
        return "figures/cell_type_by_cluster.png"

    def plot_umap_by_cluster(self,tenx, pcs=50):
        if not os.path.exists("figures/umap_by_cluster.png"):
            self.workflow.transform (
                name = "umap_by_cluster",
                func = plotting.umap_by_cluster,
                args = (
                    pypeliner.managed.InputFile(self.sce),
                    tenx,
                    self.prefix,
                    pcs
                )
            )
        return "figures/umap_by_cluster.png"


    def run_scviz(self, perplexity, components):
        if components == None:
            components = "All"
        main_fig = "perplexity_{}_regularizer_0.001_batch_size_512_learning_rate_0.01_latent_dimension_2_activation_ELU_seed_1_iter_3000.png".format(perplexity)
        scviz_dir = os.path.join(self.output, "{}_{}".format(perplexity,components))
        scviz_file = os.path.join(scviz_dir,main_fig)
        if not os.path.exists(scviz_file):
            self.workflow.transform (
                name = "{}_scviz_{}_{}".format(self.prefix, perplexity, components),
                func = SCViz.train,
                args = (
                    self.output,
                    perplexity,
                    components,
                    pypeliner.managed.InputFile(self.sce)
                )
            )

    def map_scviz(self, scviz_embedding):
        self.workflow.transform (
            name = "{}_scviz".format(self.prefix),
            func = SCViz.map,
            args = (
                pypeliner.managed.InputFile(self.sce),
                scviz_embedding,
                self.output
            )
        )

    def get_workflow(self):
        return self.workflow
