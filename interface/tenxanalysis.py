import os
import collections
import glob
import scanpy.api as sc
from utils import config
import pandas

class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_gene_bc_matrices')
        self.clustering = os.path.join(self.path, 'analysis/clustering')
        self.matrix = os.path.join(self.path, "")
        self.projection = os.path.join(self.path, 'analysis/pca/10_components/projection.csv')
        self.cellranger_tsne = os.path.join(self.path, 'analysis/tsne/2_components/projection.csv')
        self.path = None
        #self.create_scanpy_adata()

    def filtered_matrices(self):
        matrices = os.path.join(self.filtered_gene_bc_matrices, config.build + "/")
        return matrices

    def raw_matrices(self):
        matrices = os.path.join(self.raw_gene_bc_matrices, config.build + "/")
        return matrices

    def summary(self):
        web_summary = os.path.join(self,path, "web_summary.html")
        assert os.path.exists(web_summary), "No summary report found!"
        return web_summary

    def tsne_cellranger_projection(self):
        rows = open(self.cellranger_tsne, "r").read().splitlines()
        rows.pop(0)
        projection = dict()
        for row in rows:
            row = row.split(",")
            projection[row[0]] = row[1:]
        return projection

    def create_scanpy_adata(self):
        print(self.filtered_matrices(), "PATH")
        self.adata = sc.read_10x_mtx(self.filtered_matrices(), var_names='gene_symbols', cache=True)
        self.adata.var_names_make_unique()
        self.adata.barcodes = pandas.read_csv(os.path.join(self.filtered_matrices(),'barcodes.tsv'), header=None)[0]
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata)
        sc.tl.umap(self.adata)
        sc.tl.tsne(self.adata)

    def umap(self, min_dist=None):
        if min_dist is not None:
            sc.tl.umap(self.adata, min_dist=min_dist)
        return dict(zip(list(self.adata.obs.index), list(self.adata.obsm["X_umap"])))


    def tsne(self, perplexity=None):
        """
        Return a barcode to tsne embedding
        """
        if perplexity is not None:
            sc.tl.tsne(self.adata, perplexity=perplexity)
        return dict(zip(list(self.adata.obs.index), list(self.adata.obsm["X_tsne"])))

    def clustering_labels(self, k=10):
        labels = collections.defaultdict(dict)
        print(os.path.join(self.clustering, "*/clusters.csv"))
        cluster_files = glob.glob(os.path.join(self.clustering, "*/clusters.csv"))
        print("Clustering files", cluster_files)
        for cluster_file in cluster_files:
            num_clusters = os.path.splitext(cluster_file.split("/")[-2])[0]
            print(num_clusters)
            print(cluster_file)
            cells = open(cluster_file,"r").read().splitlines()
            for cell in cells[1:]:
                barcode, cluster = cell.split(",")
                labels[num_clusters][barcode] = cluster
        # print("Labels1", labels)
        # if graphclust:
        #     exit(0)
        #     labels = labels["graphclust"]
        # elif k is not none:
        #     assert "kmeans_{}_clusters".format(k) in labels.keys(), "k not found"
        print(labels.keys())
        labels = labels["kmeans_{}_clusters".format(k)]
        print("Labels",labels)
        return labels

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

    # def __eq__(self, other):
    #     return self.path == other.path
