import os

class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_gene_bc_matrices')

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


    # def __eq__(self, other):
    #     return self.path == other.path
