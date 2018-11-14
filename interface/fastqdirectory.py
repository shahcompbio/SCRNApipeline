import random
import os
import glob
import collections
import shutil

class SampleSheet(object):

    def __init__(self, filename = None, delim=","):
        if filename:
            attributes = collections.defaultdict(list)
            rows = open(filename,"r").read().splitlines()
            attribute_names = rows.pop(0).split(delim)
            for row in rows:
                for attr, value in zip(attribute_names, row.split(",")):
                    attributes[attr].append(value)
            for attr, value in attributes.items():
                print(attr.lower(), value)
                setattr(self, attr.lower(), value)


    def __add__(self, other):
        for attr, value in other.__dict__.items():
            if attr not in self.__dict__:
                setattr(self, attr.lower(), value)
            else:
                self.__dict__[attr] += value
        return self


class FastQDirectory(object):

    def __init__(self, directory, prefix, output):
        self.path = directory
        self.id = "run_{}_lsf".format(prefix)
        self.samples = self.get_samples(directory)
        self.output = output
        self.results = os.path.join(output, "{}/outs/".format(self.id))
        self.completed = False


    def get_samples(self, directory):
        sheets = glob.glob(os.path.join(directory, "*.csv"))
        sheets += glob.glob(os.path.join(directory, "*.tsv"))
        sample_sheet = SampleSheet()
        for sheet in sheets:
            sample_sheet += SampleSheet(filename=sheet)
        return sample_sheet

    def get_fastqs(self):
        return list(glob.glob(os.path.join(self.path,"*.fastq*")))

    def has_qc(self):
        dir = os.path.join(self.output,"fastqc")
        _fastqs = self.get_fastqs()
        not_complete = []
        for fastq in _fastqs:
            sample = os.path.splitext(os.path.splitext(os.path.split(fastq)[1])[0])[0]
            html = os.path.join(dir,"{}_fastqc.html".format(sample))
            if not os.path.exists(html):
                return False
        return True

    def qc_reports(self):
        dir = os.path.join(self.output,"fastqc")
        _fastqs = self.get_fastqs()
        for fastq in _fastqs:
            sample = os.path.splitext(os.path.splitext(os.path.split(fastq)[1])[0])[0]
            yield sample, os.path.join(dir,"{}_fastqc.html".format(sample))

    def __eq__(self, other):
        return self.path == self.other

    def check_status(self):
        filtered_matrices = os.path.join(self.results, "filtered_gene_bc_matrices_h5.h5")
        print(filtered_matrices)
        self.completed = os.path.exists(filtered_matrices)
        return self.completed
