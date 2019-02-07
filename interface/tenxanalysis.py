import os
import collections
import glob
import scanpy.api as sc
from sklearn import cluster
import pandas
import subprocess
from utils.config import Configuration
import matplotlib.pyplot as plt
import numpy
import shutil
import gzip

from interface.singlecellexperiment import SingleCellExperiment

config = Configuration()

class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_feature_bc_matrix')
        if not os.path.exists(self.raw_gene_bc_matrices):
            try:
                print("Raw gene matrices not found.  Looking for aggregated suffix...")
                self.raw_gene_bc_matrices = self.raw_gene_bc_matrices + "_mex"
                assert os.path.exists(self.raw_gene_bc_matrices)
                print("Aggregated suffix found.")
            except AssertionError as e:
                print("No raw gene matrices found. " + self.raw_gene_bc_matrices)
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_feature_bc_matrix')
        if not os.path.exists(self.filtered_gene_bc_matrices):
            try:
                print("Filtered gene matrices not found.  Looking for aggregated suffix...")
                self.filtered_gene_bc_matrices = self.filtered_gene_bc_matrices + "_mex"
                assert os.path.exists(self.filtered_gene_bc_matrices)
                print("Aggregated suffix found.")
            except AssertionError as e:
                print("No filtered gene matrices found." + self.filtered_gene_bc_matrices)
        self.clustering = os.path.join(self.path, 'analysis/clustering')
        self.matrix = os.path.join(self.path, "")
        self.projection = os.path.join(self.path, 'analysis/pca/10_components/projection.csv')
        self.cellranger_tsne = os.path.join(self.path, 'analysis/tsne/2_components/projection.csv')
        self.top_level = "/".join(self.path.split("/")[:-3])
        gzipped  = glob.glob(self.filtered_gene_bc_matrices + "/*.gz")
        gzipped += glob.glob(self.raw_gene_bc_matrices + "/*.gz")
        for gz in gzipped:
            if not os.path.exists(gz.strip(".gz")):
                with gzip.open(gz,"rb") as f:
                    content = f.read().decode("utf-8")
                with open(gz.strip(".gz"),"w") as f:
                    f.write(content)
                print("Decompressed {}".format(gz))

    def filtered_matrices(self):
        matrices = os.path.join(self.filtered_gene_bc_matrices, config.build + "/")
        matrices = self.filtered_gene_bc_matrices
        return matrices

    def raw_matrices(self):
        matrices = os.path.join(self.raw_gene_bc_matrices, config.build + "/")
        matrices = self.raw_gene_bc_matrices
        return matrices

    def filtered_barcodes(self):
        barcode_file = os.path.join(self.filtered_matrices(),"barcodes.tsv")
        return open(barcode_file,"r").read().splitlines()

    def raw_barcodes(self):
        barcode_file = os.path.join(self.raw_matrices(),"barcodes.tsv")
        return open(barcode_file,"r").read().splitlines()

    def raw_genes(self):
        genes_file = os.path.join(self.raw_matrices(),"genes.tsv")
        return dict([row.split() for row in open(genes_file,"r").read().splitlines()])

    def filtered_genes(self, as_list=False):
        genes_file = os.path.join(self.filtered_matrices(),"genes.tsv")
        if as_list:
            return [row.split()[0] for row in open(genes_file,"r").read().splitlines()]
        return dict([row.split() for row in open(genes_file,"r").read().splitlines()])

    def summary(self):
        web_summary = os.path.join(self.path, "web_summary.html")
        return web_summary

    def get_genes(self, sce):
        _transcripts = sce.rowData["Symbol"]
        try:
            adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        except Exception:
            adata = sc.read_10x_mtx(self.filtered_matrices())
        transcripts = []
        for symbol in _transcripts:
            if symbol not in adata.var.index:
                symbol = symbol.replace(".","-")
                if symbol not in adata.var.index:
                    symbol = symbol.split("-")
                    symbol = "-".join(symbol[:-1]) + ".{}".format(symbol[-1])
                    if symbol not in adata.var.index:
                        symbol = symbol.split(".")[0]
            transcripts.append(symbol)
        return transcripts

    def gene_map(self, sce, original=False):
        _transcripts = sce.rowData["Symbol"]
        try:
            adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        except Exception:
            adata = sc.read_10x_mtx(self.filtered_matrices())
        transcripts = {}
        for symbol in _transcripts:
            original = symbol
            if symbol not in adata.var.index:
                symbol = symbol.replace(".","-")
                if symbol not in adata.var.index:
                    symbol = symbol.split("-")
                    symbol = "-".join(symbol[:-1]) + ".{}".format(symbol[-1])
                    if symbol not in adata.var.index:
                        symbol = symbol.split(".")[0]
            if original:
                transcripts[original] = symbol
            else:
                transcripts[symbol] = original
        return transcripts


    def create_scanpy_adata(self, sce, fast_load=True, assay="counts", high_var = False):
        barcodes = sce.colData["Barcode"]
        _transcripts = sce.rowData["Symbol"]
        try:
            adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        except Exception:
            adata = sc.read_10x_mtx(self.filtered_matrices())
        adata.var_names_make_unique()
        transcripts = []
        for symbol in _transcripts:
            if symbol not in adata.var.index:
                symbol = symbol.replace(".","-")
                if symbol not in adata.var.index:
                    symbol = symbol.split("-")
                    symbol = "-".join(symbol[:-1]) + ".{}".format(symbol[-1])
                    if symbol not in adata.var.index:
                        symbol = symbol.split(".")[0]
            transcripts.append(symbol)
        adata.barcodes = pandas.read_csv(os.path.join(self.filtered_matrices(),'barcodes.tsv'), header=None)[0]
        adata = adata[:,transcripts]
        assert set(adata.var.index) == set(transcripts), "Issues with symbol conversion."
        adata = adata[barcodes,:]
        adata.var_names_make_unique()
        var_transcripts = sc.pp.highly_variable_genes(adata, flavor="cell_ranger", inplace=False)
        assert len(var_transcripts) == len(adata.var.index)
        if high_var:
            var_transcripts = [x[0] for x in zip(adata.var.index, var_transcripts) if x[1][0] == True]
            adata = adata[:,var_transcripts]
        if not fast_load:
            sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.tsne(adata)
        return adata

    def get_scvis_dimensions(self, embedding_file):
        if embedding_file is None:
            raise AssertionError("scvis requires embedding file.")
        rows = open(embedding_file,"r").read().splitlines()
        header = rows.pop(0)
        embedding = []
        for row in rows:
            row = list(map(float, row.split("\t")[1:]))
            embedding.append(row)
        return numpy.array(embedding).reshape(2,len(rows))

    def clusters(self, sce, rep=None, pcs=50, embedding_file=None, copy = False):
        adata = self.create_scanpy_adata(sce, fast_load=True)
        if rep == None or rep=="PCA":
            projection = sce.getReducedDims("PCA",n=pcs)
            adata.obsm["X_pca"] = projection.T
            sc.pp.neighbors(adata, use_rep="X_pca")
        elif rep=="TSNE":
            projection = sce.getReducedDims("TSNE")
            adata.obsm["X_tsne"] = projection.T
            sc.pp.neighbors(adata, use_rep="X_tsne")
        elif rep=="SCVIS":
            if embedding_file is None:
                raise AssertionError("scvis requires embedding file.")
            rows = open(embedding_file,"r").read().splitlines()
            header = rows.pop(0)
            embedding = []
            for row in rows:
                row = list(map(float, row.split("\t")[1:]))
                embedding.append(row)
            counts = numpy.array(embedding)
            adata.obsm["X_scvis"] = counts
            sc.pp.neighbors(adata,use_rep="X_scvis")
        sc.tl.leiden(adata, resolution=config.resolution)
        if not copy:
            return adata.obs["leiden"]
        else:
            return adata

    def markers(self, sce, n=100, pcs=50):
        adata = self.clusters(sce, pcs=pcs,copy=True)
        sc.tl.rank_genes_groups(adata, 'leiden', n_genes=n)
        markers = []
        for genes, pvals in zip(adata.uns["rank_genes_groups"]["names"],adata.uns["rank_genes_groups"]["pvals_adj"]):
            for gene, pval in zip(genes,pvals):
                if pval < 0.01:
                    markers.append(gene)
        return list(set(markers))

    def markers_by_clusters(self, sce, rep=None, pcs=50, embedding_file=None):
        adata = self.clusters(sce, pcs=pcs, embedding_file=embedding_file, copy=True)
        sc.tl.rank_genes_groups(adata, 'leiden')
        sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save='_{}_{}.png'.format(rep,pcs))
        return adata.uns

    def filtered_h5(self):
        return os.path.join(self.path, "filtered_gene_bc_matrices_h5.h5")

    def raw_h5(self):
        return os.path.join(self.path, "raw_gene_bc_matrices_h5.h5")

    def filtered_mtx(self, genes, barcodes):
        matrix = os.path.join(self.filtered_matrices(),"matrix.mtx")
        rows = open(matrix,"r").read().splitlines()
        print(rows.pop(0))
        print(rows.pop(0))
        lgenes, lcells, lnzero = rows.pop(0).split()
        print(lgenes, lcells)
        sparse_matrix = collections.defaultdict(dict)
        count = 0
        for row in rows:
            row, col, val = row.split()
            gene = genes[int(row)-1]
            barcode = barcodes[int(col)-1]
            sparse_matrix[gene][barcode] = int(val)
            count += 1
        return sparse_matrix

    @staticmethod
    def read_mtx_csv(path):
        data = dict()
        genes = open(path,"r").read()
        barcodes = genes.pop(0).split(",")[1:]
        print(len(barcodes))
        for gene in genes:
            umi_counts = list(map(int, gene[1:]))
            data[gene[0]] = dict(zip(barcodes,umi_counts))
        return data

    def __add__(self, other):

        self_barcodes = self.filtered_barcodes()
        other_barcodes = other.filtered_barcodes()

        self_genes = self.filtered_genes(as_list=True)
        other_genes = self.filtered_genes(as_list=True)
        self_gene_dict = self.filtered_genes()
        other_gene_dict = other.filtered_genes()

        self_sparse_mtx = self.filtered_mtx(self_genes, self_barcodes)
        other_sparse_mtx = self.filtered_mtx(other_genes, other_barcodes)

        barcodes = set(self_barcodes).union(set(other_barcodes))
        print("Reading Mtx")
        union_genes = set(self_sparse_mtx.keys()).union(set(other_sparse_mtx.keys()))
        matrix = dict()
        nonzero = 0
        print ("Processing Genes")
        from tqdm import tqdm
        for gene in tqdm(union_genes):
            self_umis = self_sparse_mtx[gene]
            other_umis = other_sparse_mtx[gene]
            row = dict()
            for barcode in barcodes:
                self_count = 0
                other_count = 0
                if barcode in self_umis: self_count = self_umis[barcode]
                if barcode in other_umis: other_count = other_umis[barcode]
                if self_count + other_count == 0: continue
                row[barcode] = self_count + other_count
                nonzero += 1
            matrix[gene] = row
        combined_mtx = os.path.join(self.path, "combined_gene_bc_matrices")
        try:
            os.makedirs(combined_mtx)
        except OSError:
            pass
        mtx = os.path.join(combined_mtx, "matrix.mtx")
        genes = os.path.join(combined_mtx, "genes.tsv")
        cells = os.path.join(combined_mtx, "barcodes.tsv")
        output = open(mtx,"w")
        output.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        output.write("{} {} {}\n".format(len(union_genes),len(barcodes),nonzero))
        all_genes = matrix.keys()
        row = 0
        col = 0
        use_genes = []
        use_barcodes = []
        for gene in tqdm(union_genes):
            for barcode in barcodes:
                try:
                    output.write("{} {} {}\n".format(row, col, matrix[gene][barcode]))
                    row += 1
                    col += 1
                except Exception as e:
                    continue
        output.close()
        output = open(genes,"w")
        for gene in use_genes:
            symbol = self_gene_dict.get(gene,"---")
            if symbol == "---": symbol = other_gene_dict.get(gene,"---")
            output.write("{} {}\n".format(gene,symbol))
        output.close()
        output = open(cells, "w")
        for barcode in use_barcodes:
            output.write(barcode +"\n")
        output.close()
        return combined_mtx
