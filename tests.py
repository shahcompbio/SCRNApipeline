import unittest
from interface.singlecellexperiment import SingleCellExperiment
from software.dropletutils import DropletUtils
from software.cellassign import CellAssign
from software.clonealign import CloneAlign
from software.cellranger import CellRanger
from utils.plotting import PlotRanks
from pstats import Stats
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
        pass

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
        #TODO performance issue on large assay
        #assert sce_from_rs4.assay("counts").shape == sce_from_rs4.assays["counts"].shape

    @unittest.skip("Skipping...")
    def test_annotation_rdata(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        sce_from_rdata.annotate()

    @unittest.skip("Skipping...")
    def test_barcode_ranks(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.barcodeRanks(assay)

    @unittest.skip("Skipping...")
    def test_call_empty_drops(self):
        rdata = os.path.join(base_dir, "tests/example_sce.RData")
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.emptyDrops(assay)

    @unittest.skip("Skipping...")
    def test_cell_calling_diagnostics(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_qc_metrics(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_gene_exp(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_normalize_cell_bias(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_mean_var_trend(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_dim_reduction_pca(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_dim_reduction_scviz(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_scviz_dim_reduction(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_marker_detect(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_cell_assign_em(self):
        ca = CellAssign()

    @unittest.skip("Skipping...")
    def test_unwrap(self):
        example_rda = os.path.join(base_dir, "tests/example_sce.rda")
        example_clonealign_fit = os.path.join(example_rda, "tests/example_clonealign_fit.rda")
        sce = SingleCellExperiment.fromRData(example_rda)
        print(type(sce))
        del sce

    @unittest.skip("Skipping...")
    def test_clone_align(self):
        example_rda = os.path.join(base_dir, "tests/example_sce.rda")
        example_clonealign_fit = os.path.join(example_rda, "tests/example_clonealign_fit.rda")
        sce = SingleCellExperiment.fromRData(example_rda)
        clonealigner = CloneAlign()
        res = clonealigner.run(sce)

    @unittest.skip("Skipping...")
    def test_cellranger_mkfastq(self):
        ranger = CellRanger()
        ranger.mkfastq()

    @unittest.skip("Skipping...")
    def test_cellranger_count(self):
        fastq_directory = "/home/ceglian/data/R262_S3"
        ranger = CellRanger()
        ranger.count(fastq_directory)


    @unittest.skip("Skipping...")
    def test_cellranger_aggregate(self):
        raise AssertionError("Not Implemented")

    @unittest.skip("Skipping...")
    def test_cellranger_dim_reduction(self):
        raise AssertionError(self)


if __name__ == '__main__':
    unittest.main()
    print("Tests Complete")
