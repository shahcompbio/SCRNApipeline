import os
from utils.plotting import *


class ScaterAnalysis(object):

    def __init__(self,output):
        self.output = output

    def barcodeRanks(self, ranks):
        assert "rank" in ranks
        assert "total" in ranks
        filename = os.path.join(self.output, "ranks.png")
        hlines = {"Inflection":ranks["inflection"],"Knee":ranks["knee"]}
        scatter(ranks["rank"], ranks["total"], "Rank", "Total UMI Count", hlines, filename)


    def emptyDrops(self, empty):
        print(empty.keys())
        print(len(empty["FDR"]))
        print(sum([empty["Total"]]))
        print(len([cell for cell in empty["FDR"] if cell <= 0.01]))

    def QCMetrics(self, metrics):
        pass
