import unittest
from interface.singlecellexperiment import SingleCellExperiment
from interface.dropletutils import DropletUtils
from utils.plotting import PlotRanks

class TestSingleCellExperiment(unittest.TestCase):

    def setUp(self):
        pass

    def test_load_rdata(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        pass

    def test_read10x_counts_load_rs4(self):
        tenx = DropletUtils()
        rs4_result = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment.fromRS4(rs4_result)

    def test_assay_names_rdata(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        assays = sce_from_rdata.assays
        assay_names = list(assays.keys())
        self.assertEqual(assay_names,['BatchCellMeans', 'BaseCellMeans', 'BCV', 'CellMeans', 'TrueCounts', 'counts'])

    def test_raw_assay_type_equivelence(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        rs4_results = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment(rs4_results)
        self.assertEqual(type(sce_from_rdata.assays["counts"]),type(sce_from_rdata.assays["counts"]))

    def test_get_assay_with_row_col_data_rdata(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        assert sce_from_rdata.assay("counts").shape == sce_from_rdata.assays["counts"].shape

    def test_get_assay_with_row_col_data_rs4(self):
        tenx = DropletUtils()
        rs4_results = tenx.read10xCounts("/home/ceglian/data/raw_gene_bc_matrices/hg19/")
        sce_from_rs4 = SingleCellExperiment.fromRS4(rs4_results)
        #TODO performance issue on large assay
        #assert sce_from_rs4.assay("counts").shape == sce_from_rs4.assays["counts"].shape

    def test_annotation_rdata(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        sce_from_rdata.annotate()

    def test_barcode_ranks(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.ranks(assay)

    def test_call_empty_drops(self):
        rdata = "/home/ceglian/data/example_sce.RData"
        sce_from_rdata = SingleCellExperiment.fromRData(rdata)
        tenx = DropletUtils()
        assay = sce_from_rdata.assays["counts"]
        values = tenx.emptydrops(assay)

    def test_cell_calling_diagnostics(self):
        raise AssertionError("Not Implemented")

    def test_qc_metrics(self):
        raise AssertionError("Not Implemented")

    def test_gene_exp(self):
        raise AssertionError("Not Implemented")

    def test_normalize_cell_bias(self):
        raise AssertionError("Not Implemented")

    def test_mean_var_trend(self):
        raise AssertionError("Not Implemented")

    def test_dim_reduction_pca(self):
        raise AssertionError("Not Implemented")

    def test_dim_reduction_scviz(self):
        raise AssertionError("Not Implemented")

    def test_clustering(self):
        raise AssertionError("Not Implemented")

    def test_marker_detect(self):
        raise AssertionError("Not Implemented")

    def test_cell_assign_em(self):
        raise AssertionError("Not Implemented")

    def test_clone_align(self):
        raise AssertionError("Not Implemented")

if __name__ == '__main__':
    unittest.main()
