import os
import collections
import glob

class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_gene_bc_matrices')
        self.clustering = os.path.join(self.path, 'analysis/clustering')
        self.projection = os.path.join(self.path, 'analysis/pca/10_components/projection.csv')

    def filtered_matrices(self, species):
        matrices = os.path.join(self.filtered_gene_bc_matrices, species)
        return matrices

    def raw_matrices(self, species):
        matrices = os.path.join(self.raw_gene_bc_matrices, species)
        return matrices

    def summary(self):
        web_summary = os.path.join(self,path, "web_summary.html")
        assert os.path.exists(web_summary), "No summary report found!"
        return web_summary

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
