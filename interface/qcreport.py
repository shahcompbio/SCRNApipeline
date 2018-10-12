import os
import glob

class QCReport(object):

    def __init__(self):
        self.fastq_htmls = None
        self.fastq_zipped = None
        self.has_fastqc = False
        self.fastq_directory = None
        self.has_cellranger = False
        self.cellranger_summary = None

    def add_fastqc(self, directory):
        html_files = os.path.join(directory,"*.html")
        zipped_directory = os.path.join(directory,"*.zip")
        self.fastqc_htmls = glob.glob(html_files)
        self.fastqc_zipped = glob.glob(zipped_directory)
        self.fastq_directory = directory
        self.has_fastqc = True

    def add_cellranger(self, tenx_analysis):
        self.cellranger_summary = tenx_analysis.summary()
        self.has_cellranger = True


    def add_scater(self, scater_analysis):
        pass

    def generate_report(self):
        pass

    def __eq__(self, other):
        if self.has_fastq != other.has_fastqc:
            return False
        if self.has_fastq and other.has_fastq:
            if self.fastq_directory != other.fastq_directory:
                return False
        if self.has_cellranger != other.has_cellranger:
            return False
        if self.has_cellranger and other.has_cellranger:
            if self.cellranger_summary != other.cellranger_summary:
                return False
        return True
