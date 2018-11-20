from interface.singlecellexperiment import SingleCellExperiment
from interface.tenxanalysis import TenxAnalysis
from interface.cells import Cell, Cells
import gc

# sce = SingleCellExperiment.fromRData("~/project_data/pre_37/sce_final.rdata")
# print(sce.colData.keys())
# # print(sce.rowData.keys())
exp = TenxAnalysis("/Users/ceglian/project_data/pre_37/run_pre_37_0_lsf/outs/")
adata = exp.create_scanpy_adata()
transcripts = list(adata.var.index)
barcodes = list(adata.obs.index)
data = adata.X
cells = Cells()

print(len(barcodes))
print(len(transcripts))
print(data.shape)
for barcode, i in zip(barcodes,range(data.shape[0])):
    expression = dict()
    row = list(data.getrow(i).todense().transpose())
    assert len(row) == len(transcripts)
    cell = Cell()
    cell.barcode = barcode
    cell.data = dict(zip(transcripts, row))
    cells.add_cell(cell)
    del row
    gc.collect()
