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

    def _generate_csv(self):
        csv = open(self.libraries_csv,"w")
        csv.write("library_id,molecule_h5\n")
        for library in self.libraries:
            lib_output = os.path.join(output,library)
            partial_path = os.path.join(lib_output,"*/outs/molecule_info.h5")
            for molecule_h5 in glob.glob(partial_path):
                if os.path.exists(molecule_h5):
                    csv.write("{prefix},{path}\n".format(prefix=library,path=molecule_h5))
                else:
                    print("Not Found", molecule_h5)
        csv.close()

    def aggregate_libraries(self, libraries):
        if len(libraries) > 0:
            self.libraries = libraries
            self.libraries_csv = os.path.join(self.output, "libraries.csv")
            self._generate_csv()
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
        self.sce_filtered = os.path.join(self.output,"sce_filtered.rdata")
        self.sce_final = os.path.join(self.output,"sce_final.rdata")
        self.cell_assign_fit = os.path.join(self.output, "cell_assign_fit.pkl")
        self.cell_assign_rdata = os.path.join(self.output, "{}_cellassign.rdata".format(self.prefix))
        self.clone_align_fit = os.path.join(self.output, "clone_align_fit.rdata")


    def set_directory(self, directory):
        if directory is not None:
            self.analysis = TenxAnalysis(directory)
            self.build_sce()

    def set_rdata(self, rdata):
        if rdata is not None:
            self.rdata = rdata

    def build_sce(self):
        self.workflow.transform (
            name = "{}_read10xcounts".format(self.prefix),
            func = TenX.read10xCountsFiltered,
            args = (
                self.analysis,
                pypeliner.managed.OutputFile(self.sce_filtered)
            )
        )
        self.rdata = self.sce_filtered

    def run_scater(self):
        self.scater_workflow = ScaterCode(self.output)
        self.rscript = self.scater_workflow.generate_script()
        if not os.path.exists(self.sce_final):
            self.workflow.commandline (
                name = "{}_scater_workflow".format(self.prefix),
                args = ("Rscript", self.rscript, pypeliner.managed.InputFile(self.sce_filtered), pypeliner.managed.OutputFile(self.sce_final)),
            )


    def run_cell_assign(self, rho_matrix, additional=None):
        self.workflow.transform (
            name = "{}_cellassign".format(self.prefix),
            func = CellAssign.run_em,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                pypeliner.managed.OutputFile(self.cell_assign_fit),
                self.prefix,
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

    def plotting_workflow(self):
        self.plot_cell_types()
        self.plot_tsne_by_cluster()
        self.plot_tsne_by_cell_type()
        self.plot_cell_type_by_cluster()
        self.plot_umap()

    def plot_cell_types(self):
        self.workflow.transform (
            name = "celltypes",
            func = plotting.celltypes,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                pypeliner.managed.InputFile(self.cell_assign_fit),
                self.prefix
            )
        )
        return "cell_types.png"

    def plot_tsne_by_cluster(self):
        self.workflow.transform (
            name = "tsne_by_cluster",
            func = plotting.tsne_by_cluster,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                pypeliner.managed.InputFile(self.cell_assign_fit),
                self.prefix
            )
        )
        return "tsne_by_cluster.png"

    def plot_tsne_by_cell_type(self):
        self.workflow.transform (
            name = "tsne_by_cell_type",
            func = plotting.tsne_by_cell_type,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                pypeliner.managed.InputFile(self.cell_assign_fit),
                self.prefix
            )
        )
        return "tsne_by_cell_type.png"

    def plot_cell_type_by_cluster(self):
        self.workflow.transform (
            name = "cell_type_by_cluster",
            func = plotting.cell_type_by_cluster,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                pypeliner.managed.InputFile(self.cell_assign_fit),
                self.prefix
            )
        )
        return "cell_type_by_cluster.png"

    def plot_umap(self):
        self.workflow.transform (
            name = "{}_umap".format(self.prefix),
            func = plotting.umap,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                self.analysis,
                pypeliner.managed.InputFile(self.cell_assign_fit),
                "umap.png"
            )
        )

    def run_scviz(self, cellassignments):
        self.workflow.transform (
            name = "{}_scviz".format(self.prefix),
            func = SCViz.train,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                self.output,
                cellassignments
            )
        )

    def map_scviz(self, scviz_embedding):
        workflow.transform (
            name = "{}_scviz".format(self.prefix),
            func = SCViz.map,
            args = (
                pypeliner.managed.InputFile(self.sce_final),
                scviz_embedding,
                self.output
            )
        )

    def get_workflow(self):
        return self.workflow
