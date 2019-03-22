import scanorama
import pickle
import scanpy.api as sc
import matplotlib.pyplot as plt
from interface.tenxanalysis import TenxAnalysis
from utils.cloud import TenxDataStorage


class Scanorama(object):

    @staticmethod
    def get_genes(samples):
        genes_by_sample = []
        for sample in samples:
            sce = sample.qcd_sce()
            genes_by_sample.append(set(sce.rowData["Symbol"]))
        return set.intersection(*genes_by_sample)

    @staticmethod
    def get_tenx(samples):
        tenxs = []
        for sample in samples:
            tenx = TenxDataStorage(sample)
            tenx.download()
            tenxs.append(TenxAnalysis(tenx.tenx_path))
        return tenxs

    @staticmethod
    def get_qcd_adata(tenx,subset):
        sce = tenx.qcd_adata(subset=subset)
        return sce

    @staticmethod
    def integrate_and_correct(samples, assay="counts"):
        tenxs = Scanorama.get_tenx(samples)
        common_genes = Scanorama.get_genes(tenxs)
        adatas = [Scanorama.get_qcd_adata(sample, common_genes) for sample in tenxs]
        integrateds, correcteds = scanorama.correct_scanpy(adatas, return_dimred=True)
        for tenx, integrated, corrected in zip(tenxs, integrateds, correcteds):
            tenx.set_corrected(corrected)
            tenx.set_integrated(integrated)


    @staticmethod
    def plot_integrated(tenxs):
        for tenx in tenxs:
            integrated = tenx.get_integrated()
            print (integrated[1].shape)

    @staticmethod
    def plot_corrected(tenx, sample="qcd", subset=None):
        adata = tenx.qcd_adata(subset=subset)
        corrected_adata = tenx.get_corrected(adata)
        batch_adata = adata.concatenate(corrected_adata,batch_key="correction")
        batch_adata = sc.tl.pca(batch_adata, copy=True)
        batch_adata = sc.pp.neighbors(batch_adata,copy=True)
        batch_adata = sc.tl.umap(batch_adata,copy=True)
        plt.figure()
        sc.pl.umap(batch_adata,color="correction",show=False)
        plt.savefig("{}_corrected.png".format(sample))


def main():
    samples = ["TENX003_SA605X3XB01966_002", "TENX003_SA605X3XB01966_001"]
    # tenxs = Scanorama.get_tenx(samples)
    # genes = Scanorama.get_genes(tenxs)
    # samples_to_tenx = zip(samples,tenxs)
    Scanorama.read_molecules("molecule_info.h5")
    # Scanorama.plot_integrated(tenxs)
    # for sample, tenx in samples_to_tenx:
    #     Scanorama.plot_corrected(tenx, sample=sample, subset=genes)


if __name__ == '__main__':
    main()
