from interface.tenxanalysis import TenxAnalysis
from utils.config import Configuration

import os
import pickle
import glob
import shutil

config = Configuration()

class Results(object):

    def __init__(self, output):
        self.plots = list()
        self.output = output
        self.report_dir = os.path.join(self.output,"{}_report/".format(config.prefix))
        self.sce_dir = os.path.join(self.output, "rdata")
        self.tables_dir = os.path.join(self.output, "tables")
        self.figures_dir = os.path.join(self.output, "figures")
        try:
            os.makedirs(self.report_dir)
        except Exception as e:
            pass
        try:
            os.makedirs(self.sce_dir)
        except Exception as e:
            pass
        try:
            os.makedirs(self.tables_dir)
        except Exception as e:
            pass
        try:
            os.makedirs(self.figures_dir)
        except Exception as e:
            pass
        self.paths = []

    def qc_reports(self):
        for html in glob.glob(os.path.join(self.output, "fastqc/*/*.html")):
            yield html

    def add_analysis(self, tenx):
        self.analysis = TenxAnalysis(tenx)
        summary = self.analysis.summary()
        dest = os.path.join(self.report_dir, "summary.html")
        self.paths.append((summary,dest))
        self.summary = "summary.html"

    def add_workflow(self, script):
        dest = os.path.join(self.report_dir, os.path.split(script)[1])
        self.script = os.path.split(dest)[1]
        self.paths.append((script,dest))

    def add_sce(self, sce):
        dest =  os.path.join(self.report_dir, os.path.split(sce)[1])
        self.filtered_sce = os.path.split(dest)[1]
        self.paths.append((sce,dest))

    def add_cellassign_pkl(self, pkl):
        dest =  os.path.join(self.report_dir, os.path.split(pkl)[1])
        self.pkl = os.path.split(dest)[1]
        self.paths.append((pkl,dest))

    def add_cellassign_raw(self, raw):
        dest =  os.path.join(self.report_dir, os.path.split(raw)[1])
        self.raw = os.path.split(dest)[1]
        self.paths.append((raw,dest))


    def add_plot(self, path, header, desc=""):
        plot = dict()
        dest =  os.path.join(self.report_dir, os.path.split(path)[1])
        self.paths.append((path, dest))

        plot["path"] = os.path.split(dest)[1]
        plot["header"] = header
        plot["desc"] = desc
        self.plots.append(plot)

    def finalize(self):
        for source, dest in self.paths:
            try:
                shutil.copyfile(source, dest)
            except Exception as e:
                continue

    def barcode_to_celltype(self):
        tsv = os.path.join(self.output,"barcode_to_celltype.tsv")
        output = open(tsv,"w")
        cell_assignments = pickle.load(open(self.pkl,"rb"))
        for barcode, ctype in zip(cell_assignments["Barcode"],cell_assignments["cell_type"]):
            output.write("{}\t{}\n".format(barcode,ctype))
        output.close()
        return tsv
