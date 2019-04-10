import random
import os
import glob
import collections
import shutil
import time
import subprocess

from utils.cloud import FastqDataStorage


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

    def __init__(self, directory, prefix, output, datapath=None):
        if datapath is not None:
            fastqdir = os.path.join(datapath, directory)
        else:
            fastqdir = directory
        if not os.path.exists(fastqdir):
            self.storage = FastqDataStorage(prefix)
            self.storage.set_data_path(datapath)
            directory = self.storage.download_fastqs()
        self.path = directory
        self.id = prefix
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

    def get_fastqs(self, index=True):
        fastqs = list(glob.glob(os.path.join(self.path,"*.fastq*")))
        if not index:
            fastqs = [fastq for fastq in fastqs if "I1" not in fastq]
        return fastqs

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
        self.completed = os.path.exists(filtered_matrices)
        if not self.completed:
            filtered_matrices = os.path.join(self.results, "filtered_feature_bc_matrix_h5.h5")
            self.completed = os.path.exists(filtered_matrices)
        return self.completed

    def concatenate(self,compressed=True):
        t0=time.time()
        concatenated = os.path.join(self.output,"{}_cat.fastq.qz".format(self.id))
        cmd = ["cat"]
        for fastq in self.get_fastqs(index=False):
            cmd.append(fastq)
        cmd.append(">")
        cmd.append(concatenated)
        print(" ".join(cmd))
        subprocess.call(cmd, shell=True)
        t1=time.time()
        print(t1-t0, "sec to combine.")
        if not compressed:
            subprocess.call(["gunzip","-f",concatenated])
        t2=time.time()
        print(t1-t2, "sec to decompress.")
        return concatenated
