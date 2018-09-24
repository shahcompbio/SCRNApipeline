import os

class TenxAnalysis(object):

    def __init__(self, directory):
        self.path = directory
        self.raw_gene_bc_matrices = os.path.join(self.path, 'raw_gene_bc_matrices')
        self.filtered_gene_bc_matrices = os.path.join(self.path, 'filtered_gene_bc_matrices')

    def filtered_gene_bc_matrices(self, species):
        matrices = os.path.join(self.filtered_gene_bc_matrices, species)
        return matrices

    def raw_gene_bc_matrices(self, species):
        matrices = os.path.join(self.raw_gene_bc_matrices, species)
        return matrices

    def __eq__(self, other):
        return self.path == self.other
