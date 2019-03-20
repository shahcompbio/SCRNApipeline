import h5py
import os
import glob
import subprocess
import numpy as np
import time
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
import sys
import collections
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
from multiprocessing import Pool
import pickle


from utils.cloud import TenxDataStorage
from interface.tenxanalysis import TenxAnalysis
from interface.fastqdirectory import FastQDirectory

class Kallisto(object):

    def __init__(self, fastqs, tenx, chem="v2"):
        self.tenx = tenx
        self.fastqs = fastqs
        self.chem = chem
        self.binary = "kallisto"
        self.index = None
        self.output = os.path.join(self.tenx.directory, "kallisto")
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.nthreads = 64
        self.tcc_output = os.path.join(self.output,"tcc")
        if not os.path.exists(self.tcc_output):
            os.makedirs(self.tcc_output)
        self.matrix_ec = os.path.join(self.tcc_output, "matrix.ec")
        self.matrix_tsv = os.path.join(self.tcc_output, "matrix.tcc")
        self.index = "/igo_large/reference/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx"
        self.transcript_to_gene = "/igo_large/reference/t2g.txt"
        self.sorted_bus_text = os.path.join(self.tcc_output, "output.tsv")
        self.matrix_dat = os.path.join(self.tcc_output,"tcc_matrix.dat")
        self.bus_output = os.path.join(self.tcc_output,"output.bus")
        self.sorted_bus = os.path.join(self.tcc_output,"sorted.bus")
        self.bus_matrix = os.path.join(self.tcc_output,"matrix.tsv")
        self.bustools = "bustools"
        self.transcript_to_ec = collections.defaultdict(set)
        self.gene_to_transcript = collections.defaultdict(set)
        self.gene_to_ec = collections.defaultdict(set)
        self.transcripts = os.path.join(self.tcc_output, "transcripts.txt")
        self.matrix = os.path.join(self.tcc_output,"design_matrix.dat")

    def run_pseudo(self):
        if not os.path.exists(self.bus_output):
            cmd = [self.binary,"bus","-i",self.index,"-o",self.tcc_output,"-t",str(self.nthreads),"-x","10x{}".format(self.chem)]
            for fastq in self.fastqs.get_fastqs():
                cmd.append(fastq)
            print (" ".join(cmd))
            subprocess.call(cmd)
        assert os.path.exists(self.bus_output)

    def run_bus(self):
        if not os.path.exists(self.bus_matrix):
            cmd = [self.bustools, "sort", self.bus_output, "-o", self.sorted_bus]
            subprocess.call(cmd)
            assert os.path.exists(self.sorted_bus)
            cmd = [self.bustools,"text","-o",self.bus_matrix,self.sorted_bus]
            subprocess.call(cmd)
            assert os.path.exists(self.bus_matrix)


    def tcc_matrix(self):
        if not os.path.exists(self.matrix_dat):
            assert os.path.exists(self.matrix_ec), "No EC file found."
            assert os.path.exists(self.matrix_tsv), "No bus text file found."
            COOinput = np.loadtxt( self.matrix_tsv, delimiter='\t', dtype=int)
            rows,cols,data = COOinput.T
            nonzero_ec = np.unique(rows)
            map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
            map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
            TCCmatrix = coo_matrix( (data.astype(int),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) )
            NUM_OF_CELLS = TCCmatrix.shape[1]
            print("NUM_OF_CELLS =", NUM_OF_CELLS)
            T = TCCmatrix.tocsr()
            pickle.dump(T,open(self.matrix_dat,"wb"))
        else:
            T = pickle.load(open(self.matrix_dat,"rb"))
        return T

    def setup_mapping(self):
        transcript_ids = dict()
        transcripts = open(self.transcripts,"r").read().splitlines()
        for i, transcript in enumerate(transcripts):
            transcript_ids[str(i)] = transcript
        ecs = open(self.matrix_ec,"r").read().splitlines()
        self.ecs = []
        for ec in ecs:
            ecid, transcripts = ec.split()
            self.ecs.append(ecid)
            transcripts = transcripts.split(",")
            for transcript in transcripts:
                self.transcript_to_ec[transcript_ids[transcript]].add(ecid)
        genes = open(self.transcript_to_gene,"r").read().splitlines()
        for gene in genes:
            t1, t2, symbol = gene.split()
            self.gene_to_transcript[symbol].add(t1)
            self.gene_to_transcript[symbol].add(t2)
        for symbol, transcripts in self.gene_to_transcript.items():
            for transcript in transcripts:
                for ec in self.transcript_to_ec[transcript]:
                    self.gene_to_ec[symbol].add(ec)
        self.ecs = set(self.ecs)


    def design_matrix(self):
        import tqdm
        assert os.path.exists(self.bus_matrix)
        self.setup_mapping()
        print("Setup complete")
        barcodes = list()
        tccs = list()
        counts = list()
        tcc_by_cell = open(self.bus_matrix,"r")
        valid_barcodes = set(self.tenx.filtered_barcodes())
        for tcc in tcc_by_cell:
            barcode, umi, tccid, count = tcc.split("\t")
            if barcode+"-1" in valid_barcodes:
                barcodes.append(barcode)
                counts.append(count)
                tccs.append(tccid)

        ec_counts_by_cell = collections.defaultdict(dict)
        for ec, cell, count in zip(tccs,barcodes,counts):
            ec_counts_by_cell[ec][cell] = int(count)

        design_matrix = collections.defaultdict(lambda : collections.defaultdict(dict))
        for gene, ecs in tqdm.tqdm(self.gene_to_ec.items()):
            for ec in ecs:
                design_matrix[gene][ec] = ec_counts_by_cell[ec]
        return design_matrix


def main():
    sample = "patient2"

    tenx = TenxDataStorage(sample,version="v2")
    tenx.download()
    tenx_analysis = TenxAnalysis(tenx.tenx_path)
    tenx_analysis.load()
    output = "/igo_large/scratch/test_kallisto"
    fastq_directory = FastQDirectory("/igo_large/scratch/allen/bams/xfastqs2/McGilvery_Sonya__TLH_MissingLibrary_1_CB8R9ANXX/", sample, output)

    krunner = Kallisto(fastq_directory, tenx_analysis)
    krunner.de()

if __name__ == '__main__':
    main()
