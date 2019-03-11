import h5py
import os
import glob
import subprocess
import numpy as np
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
import sys
import collections
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import *
from scipy.stats import entropy

from utils.cloud import TenxDataStorage
from interface.tenxanalysis import TenxAnalysis

"""
try:
    os.system('python get_cell_barcodes.py '+json_path)
    if not os.path.isfile(str(parameter["SAVE_DIR"])+"barcodes.dat"):
        print("ERROR:"+str(parameter["SAVE_DIR"])+"barcodes.dat not found"); exit(1)
    os.system('python error_correct_and_split.py '+json_path)
    if not os.path.isfile(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt"):
        print("ERROR:"+str(parameter["OUTPUT_DIR"])+"umi_read_list.txt not found"); exit(1)
    os.system('python compute_TCCs.py '+json_path)
    if not os.path.isfile(parameter["kallisto"]["TCC_output"]+"matrix.tsv"):
        print("ERROR:"+parameter["kallisto"]["TCC_output"]+"matrix.tsv not found"); exit(1)
os.system('python prep_TCC_matrix.py '+json_path)


barcodes:
with open(save_dir+"barcodes.dat", 'wb') as f:
    pickle.dump(barcodes,f)
with open(save_dir+"codewords.dat", 'wb') as f:
    pickle.dump(codewords,f)
with open(save_dir+"brc_idx_to_correct.dat", 'wb') as f:
pickle.dump(brc_idx_to_correct,f)


error_correct_and_split
with open(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt", 'wb') as f:
f.write(out_data)

"""


class Kallisto(object):

    def __init__(self, sample, chem="v2"):
        self.chem = chem
        tenx = TenxDataStorage(sample,version=self.chem)
        tenx.download()
        self.tenx = TenxAnalysis(tenx.tenx_path)
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
        self.matrix_tsv = os.path.join(self.tcc_output, "final_matrix.tsv")
        self.index = "/igo_large/reference/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx"
        self.transcript_to_gene = "/igo_large/reference/t2g.txt"
        self.read_molecules()
        self.fastq_dir = "/igo_large/scratch/allen/bams/xfastqs2/McGilvery_Sonya__TLH_MissingLibrary_1_CB8R9ANXX/*.fastq.gz"
        self.fastqs = glob.glob(self.fastq_dir)
        self.sorted_bus_text = os.path.join(self.tcc_output, "output.tsv")
        self.matrix_dat = os.path.join(self.tcc_output,"tcc_matrix.dat")
        self.l1_distances = os.path.join(self.tcc_output, "l1_distances.dat")
        self.nonzero_ec = os.path.join(self.tcc_output, "nonzero_ec.dat")

    def run(self):
        #cmd = [self.binary,"pseudo","-i",self.index,"-o",self.tcc_output,"-t",str(self.nthreads),"-x","10x{}".format(self.chem)]
        cmd = [self.binary,"pseudo","-i",self.index,"-o",self.tcc_output,"-t",str(self.nthreads)]
        for fastq in self.fastqs:
            cmd.append(fastq)
        print (" ".join(cmd))
        subprocess.call(cmd)

    def read_molecules(self):
        assert os.path.exists(self.tenx.molecules_h5()), self.tenx.molecules_h5()
        self.molecules = h5py.File(self.tenx.molecules_h5(), "r")
        self.valid_keys = list(self.molecules.keys())
        #barcodes
        #reads per barcode

    def prep_matrix(self):
        assert os.path.exists(self.matrix_ec), "No EC file found."
        assert os.path.exists(self.matrix_tsv), "No EC file found."
        COOinput = np.loadtxt( self.matrix_tsv, delimiter='\t' , dtype=int)
        rows,cols,data = COOinput.T
        nonzero_ec = np.unique(rows)
        map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
        map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
        TCCmatrix = coo_matrix( (data.astype(float),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) )

    def de(self):
        pass

    def generate_matrix_from_bus(self):
        gene_min = 200
        gene_max = 10000
        tr2g = {}
        trlist = []
        with open(self.transcript_to_gene) as f:
            for line in f:
                l = line.split()
                tr2g[l[0]] = l[2]
                trlist.append(l[0])

        genes = list(set(tr2g[t] for t in tr2g))

        # load equivalence classes
        ecs = {}
        with open(self.matrix_ec) as f:
            for line in f:
                l = line.split()
                ec = int(l[0])
                trs = [int(x) for x in l[1].split(',')]
                ecs[ec] = trs

        def ec2g(ec):
            if ec in ecs:
                return list(set(tr2g[trlist[t]] for t in ecs[ec]))
            else:
                return []

        cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
        pbar=None
        pumi=None
        with open(self.sorted_bus_text) as f:
            gs = set()
            for line in f:
                l = line.split()
                try:
                    barcode,umi,ec,count = line.split()
                except Exception as e:
                    continue
                ec = int(ec)

                if barcode == pbar:
                    # same barcode
                    if umi == pumi:
                        # same UMI, let's update with intersection of genelist
                        gl = ec2g(ec)
                        gs.intersection_update(gl)
                    else:
                        # new UMI, process the previous gene set
                        for g in gs:
                            cell_gene[barcode][g] += 1.0/len(gs)
                        # record new umi, reset gene set
                        pumi = umi
                        gs = set(ec2g(ec))
                else:
                    # work with previous gene list
                    for g in gs:
                        cell_gene[pbar][g] += 1.0/len(gs)

                    if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                        del cell_gene[pbar]

                    pbar = barcode
                    pumi = umi

                    gs = set(ec2g(ec))

            for g in gs:
                cell_gene[pbar][g] += 1.0/len(gs)

            if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                del cell_gene[pbar]

        barcode_hist = collections.defaultdict(int)
        for barcode in cell_gene:
            cg = cell_gene[barcode]
            s = len([cg[g] for g in cg])
            barcode_hist[barcode] += s

        #Output a gene count histogram
        bcv = [x for b,x in barcode_hist.items() if x > gene_min and x < gene_max]
        plt.switch_backend('agg')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(bcv,bins=100)
        ax.set_title("Histogram")
        plt.xlabel("number of genes detected")
        plt.ylabel("number of barcodes")
        fig.savefig('gene_hist.png')

        outfile = self.matrix_tsv

        gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
        barcodes_to_use = [b for b,x in barcode_hist.items() if x > gene_min and x < gene_max]

        num_entries = 0
        for barcode in barcodes_to_use:
            num_entries += len([x for x in cell_gene[barcode].values() if x>0])

        with open(outfile, 'w') as of:
            of.write('%%MatrixMarket matrix coordinate real general\n%\n')
            #number of genes
            of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), round(num_entries)))
            bcid = 0
            for barcode in barcodes_to_use:
                bcid += 1
                cg = cell_gene[barcode]
                gl = [(gene_to_id[g],cg[g]) for g in cg if cg[g] > 0]
                gl.sort()
                for x in gl:
                    of.write("%d %d %f\n"%(x[0],bcid,x[1]))

        gene_names = {}
        with open(self.transcript_to_gene) as f:
            f.readline()
            for line in f:
                t,g,gn = line.split()
                gene_names[g] = gn

        id_to_genes = dict((i,g) for (g,i) in gene_to_id.items())
        gl = []
        for i in range(1,len(genes)+1):
            g = id_to_genes[i]
            gid = g
        #    gid = g[:g.find('.')]
            if gid in gene_names:
                gn = gene_names[gid]
            else:
                gn = ''
            gl.append((g,gn))

        with open('./genes.tsv','w') as of:
            for g,gn in gl:
                of.write("%s\t%s\n"%(g,gn))

        with open('./barcodes.tsv','w') as of:
            of.write('\n'.join(x + '' for x in barcodes_to_use))
            of.write('\n')

    def prep_tcc_matrix(self):
        #matrix.ec file
        ecfile_dir = self.matrix_ec
        tsvfile_dir = self.matrix_tsv

        print("Loading TCCs..")

        COOinput = np.loadtxt( self.matrix_tsv, dtype=float)
        rows,cols,data = COOinput.T
        nonzero_ec = np.unique(rows)
        map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
        map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
        TCCmatrix   = coo_matrix( (data.astype(float),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) )

        NUM_OF_CELLS = TCCmatrix.shape[1]
        print("NUM_OF_CELLS =", NUM_OF_CELLS)

        T = TCCmatrix.tocsr()
        T_norm = normalize(T, norm='l1', axis=0)
        T_normT = T_norm.transpose()
        del TCCmatrix;
        _ = gc.collect()


        print("Loading Barcodes...")
        t0 = time.time()
        with open(save_dir+"barcodes.dat", 'rb') as f:
            barcodes=pickle.load(f)
        with open(save_dir+"codewords.dat", 'rb') as f:
            codewords=pickle.load(f)
        with open(save_dir+"brc_idx_to_correct.dat", 'rb') as f:
        brc_idx_to_correct= pickle.load(f)


        # Pairwise_distances

        def L1_distance(p,q):
            return cityblock(p,q).sum()

        # def jensen_shannon(p, q):
        #     m=0.5*p+0.5*q
        #     p = np.transpose(p[p > 0])
        #     q = np.transpose(q[q > 0])
        #     m = np.transpose(m[m > 0])
        #     return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))

        num_of_threads = self.nthreads
        print("Calculating pairwise L1 distances... ( num_threads =",num_of_threads,")")

        # D_js = pairwise_distances(T_normT,metric=jensen_shannon,n_jobs=num_of_threads)
        D_l1 = pairwise_distances(T_normT,metric=L1_distance,n_jobs=self.nthreads)

        print("writing data...")

        #Save data
        import pickle

        with open(self.matrix_dat, 'wb') as f:
            pickle.dump(T,f)
        with open(self.l1_distances, 'wb') as f:
            pickle.dump(D_l1,f)
        with open(self.nonzero_ec, 'wb') as f:
            pickle.dump(nonzero_ec,f)

        print("DONE.")

    def generate_scrna():
        pass


def main():
    samples = ["patient2"]
    krunner = Kallisto(samples[0])
    krunner.prep_tcc_matrix()

if __name__ == '__main__':
    main()
