import os

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

    def clustering_labels(self, graphclust=False, k=10):
        labels = collections.defaultdict(dict)
        cluster_files = glob.glob(os.path.join(self.clustering, "*/clusters.tsv"))
        for cluster_file in cluster_files:
            num_clusters = cluster_file.split("/")[-1].strip(".csv")
            cells = open(cluster_file,"r").read().splitlines()
            for cell in cells[1:]:
                barcode, cluster = cell.split(",")
                labels[num_clusters][barcode] = cluster
        if graph_clust:
            labels = labels["graphclust"]
        elif k is not none:
            labels = labels["kmeans_{}_clusters".format(k)]
        return labels

    def map_projection(self, cellassignments, output):
        output = open(output,"w")
        rows = open(self.projection,"r").read().splitlines()
        header = rows.pop(0).split(",")
        header[0] = "Cell_Type"
        output.write("\t".join(header)+"\n")
        for row in rows:
            row = row.split(",")
            row[0] = cellassignments[row[0]]
            output.write("\t".join(row)+"\n")
        output.close()
        
    # def __eq__(self, other):
    #     return self.path == other.path
