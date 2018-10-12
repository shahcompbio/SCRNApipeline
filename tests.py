import unittest
from rpy2.robjects import r, pandas2ri
from interface.singlecellexperiment import SingleCellExperiment
from interface.genemarkermatrix import GeneMarkerMatrix
from interface.tenxanalysis import TenxAnalysis
from interface.fastqdirectory import FastQDirectory
from interface.qcreport import QCReport
from software.dropletutils import DropletUtils
from software.cellassign import CellAssign
from software.clonealign import CloneAlign
from software.cellranger import CellRanger
from software.fastqc import FastQC
from software.tenx import TenX
from software.scviz import SCViz
from pstats import Stats
import numpy
import cProfile
import os
import sys
import glob

base_dir = os.path.dirname(os.path.abspath(__file__))

class TestSingleCellExperiment(unittest.TestCase):

    def setUp(self):
        self.profile = cProfile.Profile()
        self.profile.enable()

    def tearDown(self):
        p = Stats(self.profile)
        p.strip_dirs()
        p.sort_stats('cumtime')
        p.print_stats(12)


    @unittest.skip("Skipping...")
    def test_load_rdata(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)

    @unittest.skip("Skipping...")
    def test_read10x_counts_load_rs4(self):
        tenx = DropletUtils()
        rs4_result = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment.fromRS4(rs4_result)

    @unittest.skip("Skipping...")
    def test_assay_names_rdata(self):
        expected_assays = ['BatchCellMeans', 'BaseCellMeans', 'BCV', 'CellMeans', 'TrueCounts', 'counts']
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        assays = sce_from_rdata.assays
        assay_names = list(assays.keys())
        for assay in assay_names:
            self.assertTrue(assay in expected_assays)

    @unittest.skip("Skipping...")
    def test_raw_assay_type_equivelence(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        rs4_results = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment(rs4_results)
        self.assertEqual(type(sce_from_rdata.assays["counts"]),type(sce_from_rdata.assays["counts"]))

    @unittest.skip("Skipping...")
    def test_get_assay_with_row_col_data_rdata(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        assert sce_from_rdata.assay("counts").shape == sce_from_rdata.assays["counts"].shape

    @unittest.skip("Skipping...")
    def test_get_assay_with_row_col_data_rs4(self):
        tenx = DropletUtils()
        rs4_results = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment.fromRS4(rs4_results)
        assert sce_from_rs4.assay("counts").shape == sce_from_rs4.assays["counts"].shape

    @unittest.skip("Skipping...")
    def test_barcode_ranks(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = TenX()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.barcodeRanks(assay)
        print(values.keys())

    @unittest.skip("Skipping...")
    def test_call_empty_drops(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = TenX()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.emptyDrops(assay)
        print(values.keys())

    @unittest.skip("Skipping...")
    def test_expression_normalization(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = TenX()
        assay = sce_from_rdata.assays["counts"]
        cpm = tenx.calculateCPM(assay)
        tpm = tenx.calculateTPM(assay)
        fpkm = tenx.calculateFPKM(assay)
        assert cpm.shape == assay.shape
        assert tpm.shape == assay.shape
        assert fpkm.shape == assay.shape


    @unittest.skip("Skipping...")
    def test_load_qc(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce = SingleCellExperiment.fromRData(rdata)
        tenx = TenX()
        sce = tenx.calculateQCMetrics(sce)
        assert len(sce.rowData["n_cells_by_counts"]) == len(sce.rownames)

    def test_tenx_full_analysis(self):
        tenx = TenxAnalysis("tests")
        print("Reading Counts")
        sce = TenX.read10xCounts(tenx)
        # rdata = os.path.join(base_dir, "tests/example_sce.RData")
        # sce = SingleCellExperiment.fromRData(rdata)
        tenx = TenX()
        print("Generating Scater Analysis")
        scater_analysis = tenx.analysis(sce)


    @unittest.skip("Skipping...")
    def test_fastqc(self):
        output = "./tests/fastqc_test/"
        fastq_path = "/Users/ceglian/input_data/MICHELLE_0065_AHGNCGDMXX/Project_06000_EJ/Sample_cDNA_Pre_IGO_06000_EJ_1/"
        fqobj = FastQDirectory(fastq_path)
        fastqc = FastQC()
        fastqc.run(fqobj, output)

    @unittest.skip("Skipping...")
    def test_load_qc_report_obj(self):
        qcreport = QCReport()
        qcreport.add_fastqc(fastqc_output)
        qcreport.add_cellranger(tenx_analysis)
        qcreport.add_scater(scater_analysis)
        #load fastqc reports
        #load cellranger qc
        #load scater qc / normalization
        pass

    @unittest.skip("Skipping...")
    def test_load_analysis_obj(self):
        #load tenx analysis
        #load scvis embedding
        #load cellassign fit
        #load clonealign fit
        pass

    @unittest.skip("Skipping...")
    def test_generate_report(self):
        #single report
        #### qc report
        #### downstream report
        pass

    @unittest.skip("Skipping...")
    def test_installation(self):
        #see if we can load all libraries
        pass

    @unittest.skip("Skipping...")
    def test_configuration(self):
        #test find all relevent files
        #rho matrix
        #reference
        #build
        #all software (cellranger, rpackages, scvis, fastqc)
        pass

    @unittest.skip("Skipping...")
    def test_scviz_train(self):
        scviz_runner = SCViz()
        matrix = "./tests/bipolar_pca100.tsv"
        outdir = "./tests/scvis_output/"
        labels = "./tests/bipolar_label.tsv"
        scviz_runner.train(matrix, outdir, labels)

    @unittest.skip("Skipping...")
    def test_scvis_map(self):
        scviz_runner = SCViz()
        matrix = "./tests/retina_pca100_bipolar.tsv"
        outdir = "./tests/scvis_output_bipolar/"
        embedding = "./tests/scvis_output/model/perplexity_10_regularizer_0.001_batch_size_512_learning_rate_0.01_latent_dimension_2_activation_ELU_seed_1_iter_5400.ckpt"
        scviz_runner.map(matrix, outdir, embedding)


    @unittest.skip("Skipping...")
    def test_cell_assign_em(self):
        example_rda = os.path.join(base_dir, "tests/cell_assign_test.RData")
        sce = SingleCellExperiment.fromRData(example_rda)
        cellassigner = CellAssign()
        rho_matrix = dict()
        rho_matrix["Group1"] = ["Gene161", "Gene447", "Gene609", "Gene754", "Gene860", "Gene929", "Gene979"]
        rho_matrix["Group2"] = ["Gene161", "Gene447", "Gene609", "Gene754", "Gene860", "Gene929", "Gene979","Gene101","Gene212","Gene400"]
        rho = GeneMarkerMatrix(rho_matrix)
        res = cellassigner.run_em(sce, rho)

    @unittest.skip("Skipping...")
    def test_cell_assign_pkl(self):
        import pickle
        import collections
        tenx = TenxAnalysis("tests/pre_igo")
        sce = TenX.read10xCounts(tenx)
        handle = open("tests/rho_up.pkl","rb")
        rho_matrix = pickle.load(handle)
        handle.close()
        rho = GeneMarkerMatrix(rho_matrix)
        cellassigner = CellAssign()
        res = cellassigner.run_em(sce, rho)

    @unittest.skip("Skipping...")
    def test_symbol_retrieve(self):
        tenx = TenxAnalysis("tests/pre_igo")
        sce = TenX.read10xCounts(tenx)
        print(sce.rowData.keys())
        example_rda = os.path.join(base_dir, "tests/example_sce.rda")
        sce = SingleCellExperiment.fromRData(example_rda)
        print(sce.rowData.keys())
        tenx = DropletUtils()
        rs4_result = tenx.read10xCounts("tests/hg19/")
        sce = SingleCellExperiment.fromRS4(rs4_result)
        print(sce.rowData.keys())
        example_rda = os.path.join(base_dir, "tests/example_copy_number.rda")
        sce = SingleCellExperiment.fromRData(example_rda)
        print(sce.rowData.keys())
        print(sce.rownames)
        print(sce.colnames)

    @unittest.skip("Skipping...")
    def test_clone_align(self):
        example_rda = os.path.join(base_dir, "tests/example_sce.rda")
        sce = SingleCellExperiment.fromRData(example_rda)
        rowdata = sce.rowData
        cnv_data = []
        for column in ["A","B","C"]:
            column = rowdata[column]
            cnv_data.append(column)
        cnv_data = numpy.transpose(numpy.array(cnv_data))
        print(cnv_data.shape)
        clonealigner = CloneAlign()
        result = clonealigner.run(sce,cnv_data)
        assert len(result["clone"]) == 200


if __name__ == '__main__':
    unittest.main()
    print("Tests Complete")
