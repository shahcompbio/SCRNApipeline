import os
import collections
import glob
import scanpy.api as sc
from utils import config
import pandas
import subprocess
from utils import config
import numpy

from interface.singlecellexperiment import SingleCellExperiment


class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        if not os.path.exists(self.raw_gene_bc_matrices):
            try:
                print("Raw gene matrices not found.  Looking for aggregated suffix...")
                self.raw_gene_bc_matrices = self.raw_gene_bc_matrices + "_mex"
                assert os.path.exists(self.raw_gene_bc_matrices)
                print("Aggregated suffix found.")
            except AssertionError as e:
                print("No raw gene matrices found. " + self.raw_gene_bc_matrices)
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_gene_bc_matrices')
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

    def filtered_matrices(self):
        matrices = os.path.join(self.filtered_gene_bc_matrices, config.build + "/")
        return matrices

    def raw_matrices(self):
        matrices = os.path.join(self.raw_gene_bc_matrices, config.build + "/")
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
        #assert os.path.exists(web_summary), "No summary report found!"
        return web_summary

    def tsne_cellranger_projection(self):
        rows = open(self.cellranger_tsne, "r").read().splitlines()
        rows.pop(0)
        projection = dict()
        for row in rows:
            row = row.split(",")
            projection[row[0]] = row[1:]
        return projection

    def create_scanpy_adata(self, fast_load=True):
        print(self.filtered_matrices(), "PATH")
        self.adata = sc.read_10x_h5(self.filtered_h5(), genome=config.build)
        self.adata.var_names_make_unique()
        self.adata.barcodes = pandas.read_csv(os.path.join(self.filtered_matrices(),'barcodes.tsv'), header=None)[0]
        if not fast_load:
            sc.tl.pca(self.adata)
            sc.pp.neighbors(self.adata)
            sc.tl.umap(self.adata)
            sc.tl.tsne(self.adata)
        return self.adata

    def clusters(self, sce, rep=None, pcs=2, embedding_file=None):
        adata = self.create_scanpy_adata(fast_load=False)
        adata.var_names_make_unique()
        valid_transcripts = sce.rowData["Symbol"]
        valid_barcodes = sce.colData["Barcode"]
        adata = adata[valid_barcodes,:]
        adata = adata[:,valid_transcripts]
        print("Computing Neighbors...")
        if rep == None:
            sc.pp.neighbors(adata)
            sc.tl.louvain(adata,resolution=0.07)
        elif rep=="TSNE":
            embedding = sce.reducedDims["TSNE"]
            counts = []
            for i in range(0,len(embedding),pcs):
                counts.append(embedding[i:i+(pcs)])
            counts = numpy.array(counts)
            adata.obsm["TSNE"] = numpy.array(counts)
            sc.pp.neighbors(adata,n_neighbors=10,use_rep="TSNE")
            sc.tl.louvain(adata,resolution=0.05)
        elif rep=="SCVIS":
            if embedding_file is None:
                raise AssertionError("scvis requires embedding file.")
            rows = open(embedding_file,"r").read().splitlines()
            header = rows.pop(0)
            embedding = []
            for row in rows:
                row = list(map(float, row.split("\t")[1:]))
                embedding.append(row)
            embedding = numpy.array(embedding)
            print(embedding.shape)
            adata.obsm["SCVIS"] = embedding
            sc.pp.neighbors(adata,n_neighbors=10,use_rep="SCVIS")
            sc.tl.louvain(adata,resolution=0.2)
        print("Finding Clusters...")

        return adata.obs["louvain"].to_dict()

    def umap(self, min_dist=None):
        if min_dist is not None:
            sc.tl.umap(self.adata, min_dist=min_dist)
        return dict(zip(list(self.adata.obs.index), list(self.adata.obsm["X_umap"])))

    def tsne(self, perplexity=None):
        if perplexity is not None:
            sc.tl.tsne(self.adata, perplexity=perplexity)
        return dict(zip(list(self.adata.obs.index), list(self.adata.obsm["X_tsne"])))

    def clustering_labels(self, k=10):
        labels = collections.defaultdict(dict)
        cluster_files = glob.glob(os.path.join(self.clustering, "*/clusters.csv"))
        for cluster_file in cluster_files:
            num_clusters = os.path.splitext(cluster_file.split("/")[-2])[0]
            cells = open(cluster_file,"r").read().splitlines()
            for cell in cells[1:]:
                barcode, cluster = cell.split(",")
                labels[num_clusters][barcode] = cluster
        labels = labels["kmeans_{}_clusters".format(k)]
        return labels

    def get_tsne_plots_by_perplexity(self):
        if len(glob.glob(os.path.join(self.top_level, "tsne_*.png"))) == 0: exit(0)
        for plot in glob.glob(os.path.join(self.top_level, "tsne_*.png")):
            yield plot

    def map_projection(self, cellassignments, output):
        output = open(output,"w")
        rows = open(self.projection,"r").read().splitlines()
        header = rows.pop(0).split(",")
        header[0] = "Cell_Type"
        output.write("\t".join(header[1:])+"\n")
        for row in rows:
            row = row.split(",")
            try:
                cellassignments[row[0]]
                output.write("\t".join(row[1:])+"\n")
            except KeyError as e:
                continue
        output.close()

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
                    # use_genes.append(gene)
                    # use_barcodes.append(barcode)
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
