import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from interface.singlecellexperiment import SingleCellExperiment, ExpressionSet, ExpressionSetFactory

CellAssignInterface = importr("cellassign")
BiobaseInterface    = importr("Biobase")
SummarizedExperimentInterface = importr('SummarizedExperiment')
SingleCellExperimentInterface = importr('SingleCellExperiment')

class CellAssign(object):

    def __init__(self):
        pass

    def run_em(self, sce_experiment):
        # exprs = SummarizedExperimentInterface.SummarizedExperiment()
        # exprs.assays.slots["counts"] = sce_experiment.assays["counts"]
        # print(type(exprs.slots["assayData"]))
        #
        # for obj, value in robjects.__dict__.items():
        #     try:
        #         print(obj, type(value))
        #     except Exception as e:
        #         print(e)
        #         continue
        rs4 = robjects.r[r.load("tests/example_sce.RData")[0]]
        data = robjects.r["SummarizedExperiment"](rs4)
        print(SummarizedExperimentInterface.rowData(rs4))
        print(SummarizedExperimentInterface.colData(rs4))
        print(SummarizedExperimentInterface.assays(rs4))

        print("\n")
        print("Assays",data.slots["assays"].__dict__)
        print("\n")
        # for slot, value in data.slots.items():
        #     print(slot, value)

        #data.assays = SummarizedExperimentInterface.assays(data)
        #print(list(SummarizedExperimentInterface.assays(data).slotnames()))

        #data.assays = SummarizedExperimentInterface.assays(sce_experiment)
        print(CellAssignInterface.cellassign_em(data))
