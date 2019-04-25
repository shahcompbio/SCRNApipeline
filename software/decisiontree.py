from sklearn import tree
import pickle
import scanpy.api as sc

def get_rho(path):
    return None

def build(pyfit):
    celltypes = pickle.load(open(pyfit,"rb"))
    adata = sc.read_10x_mtx(tenx_path,var_names='gene_symbols')
    rho = get_rho(rho_path)
    adata =
