from interface.tenxanalysis import TenxAnalysis
import os
import pickle

class Results(object):

    def __init__(self, output, prefix):
        self.plots = list()
        self.output = output
        self.prefix = prefix

    def qc_reports(self):
        for html in glob.glob(self.output, "fastqc/*/*.html"):
            yield html

    def add_analysis(self, tenx):
        self.analysis = TenxAnalysis(tenx)

    def add_workflow(self, script):
        self.script = script

    def add_filtered_sce(self, sce):
        self.filtered_sce = sce

    def add_final_sce(self, sce):
        self.final_sce = sce

    def add_cellassign_pkl(self, pkl):
        self.pkl = pkl

    def add_cellassign_raw(self, raw):
        self.raw = raw

    def add_plot(self, path, header, desc=""):
        plot = dict()
        plot["path"] = path
        plot["header"] = header
        plot["desc"] = desc
        self.plots.append(plot)

    def barcode_to_celltype(self):
        tsv = os.path.join(self.output,"barcode_to_celltype.tsv")
        output = open(tsv,"w")
        cell_assignments = pickle.load(open(self.pkl,"rb"))
        for barcode, ctype in zip(cell_assignments["Barcode"],cell_assignments["cell_type"]):
            output.write("{}\t{}\n".format(barcode,ctype))
        output.close()
        return tsv
