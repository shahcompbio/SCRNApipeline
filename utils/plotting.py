import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import math
from sklearn import preprocessing
from utils.config import *
import pickle
from interface.genemarkermatrix import GeneMarkerMatrix
from interface.singlecellexperiment import SingleCellExperiment
from interface.tenxanalysis import TenxAnalysis
import scanpy.api as sc
import numpy
import collections
import os
import sys
from PIL import Image
import fpdf

sns.set(style="darkgrid")

def combine_figures(figures,filename):
    pdf = fpdf.FPDF()
    for figure in figures:
        pdf.add_page()
        pdf.image(figure,10,10,250,150)
    pdf.output(filename, "F")

def celltypes(rdata, cell_assign_fit, prefix, output):
    sce = SingleCellExperiment.fromRData(rdata)
    valid_barcodes = sce.colData["Barcode"]
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cell_types = list(fit["cell_type"])
    barcodes = list(fit["Barcode"])
    order = []
    invalid = 0
    for barcode, cell_type in zip(barcodes[:len(valid_barcodes)],cell_types[:len(valid_barcodes)]):
        if barcode not in valid_barcodes:
            invalid += 1
            continue
        if cell_type not in order:
            order.append(cell_type)
        if len(order) == len(set(cell_types)):
            break
    print("Dismissed {} barcodes.".format(invalid))
    f, ax = plt.subplots(figsize=(12,6))
    ax.set_title("Cell Type Assignments - {}".format(prefix))
    sns.countplot(cell_types, palette="tab10", order=list(sorted(order)))
    ax.set_xticklabels(labels=list(sorted(order)),rotation=90)
    plt.tight_layout()
    figure = os.path.join(output, "figures/cell_types.png")
    plt.savefig(figure)

def scvis_by_cell_type(rdata, cell_assign_fit, prefix, embedding_file):
    print("yas")
    fit = pickle.load(open(cell_assign_fit,"rb"))
    sce = SingleCellExperiment.fromRData(rdata)
    barcodes = sce.colData["Barcode"]
    cell_types = dict(zip(fit["Barcode"][:len(barcodes)],fit["cell_type"][:len(barcodes)]))
    rows = open(embedding_file,"r").read().splitlines()
    dims = []
    rows.pop(0)
    for row in rows:
        row = row.split("\t")
        row = list(map(float,row[1:]))
        dims.append(row)
    x = []
    y = []
    clusters = []
    for barcode, dim in zip(barcodes,dims):
        try:
            clusters.append(cell_types[barcode])
        except KeyError as e:
            clusters.append("Other")
        x.append(dim[0])
        y.append(dim[1])
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("SCVIS - Cell Type - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    print("yes")
    plt.savefig("figures/scvis_by_cell_type_{}.png".format(prefix))

def tsne_by_cell_type(rdata, cell_assign_fit, prefix):
    sce = SingleCellExperiment.fromRData(rdata)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    tsne_dims = sce.reducedDims["TSNE"]
    barcodes = sce.colData["Barcode"]
    cell_types = dict(zip(fit["Barcode"][:len(barcodes)],fit["cell_type"][:len(barcodes)]))
    tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode, cluster in cell_types.items():
        try:
            x_val = x_coded[barcode]
            y_val = y_coded[barcode]
        except Exception as e:
            continue
        try:
            clusters.append(cell_types[barcode])
        except Exception as e:
            clusters.append("Other")
        x.append(x_val)
        y.append(y_val)
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("TSNE - Cell Type - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/tsne_by_cell_type.png")

def pca_by_cell_type(rdata, cell_assign_fit, prefix):
    sce = SingleCellExperiment.fromRData(rdata)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    tsne_dims = sce.getReducedDims("PCA")
    barcodes = sce.colData["Barcode"]
    cell_types = dict(zip(fit["Barcode"][:len(barcodes)],fit["cell_type"][:len(barcodes)]))
    #tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode, cluster in cell_types.items():
        try:
            x_val = x_coded[barcode]
            y_val = y_coded[barcode]
        except Exception as e:
            continue
        try:
            clusters.append(cell_types[barcode])
        except Exception as e:
            clusters.append("Other")
        x.append(x_val)
        y.append(y_val)
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("PCA - Cell Type - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/pca_by_celltype.png")

def tsne_by_cluster(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce, pcs=pcs)
    tsne_dims = sce.reducedDims["TSNE"]
    barcodes = sce.colData["Barcode"]
    tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode, cluster in cluster_labels.items():
        clusters.append("Cluster {}".format(cluster))
        x.append(x_coded[barcode])
        y.append(y_coded[barcode])
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("TSNE - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/tsne_by_cluster.png")

def tsne_by_cluster_markers(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    markers = tenx.markers_by_clusters(sce, rep="TSNE", pcs=pcs)
    print(markers.keys())

def pca_by_cluster_markers(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    markers = tenx.markers_by_clusters(sce, rep="PCA", pcs=pcs)
    print(markers.keys())

def cluster_markers(rdata, tenx_analysis, rep, pcs, embedding_file, prefix):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    markers = tenx.markers_by_clusters(sce, rep="PCA", pcs=pcs)
    markers_by_cluster = list(zip(*markers["rank_genes_groups"]["names"]))
    for i, markers in enumerate(markers_by_cluster):
        cluster_prefix = "Cluster {} {}".format(i, prefix)
        plot_by_markers(rdata, tenx_analysis, markers, cluster_prefix, rep, pcs, embedding_file)

def umap_by_cluster_markers(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    markers = tenx.markers_by_clusters(sce, rep="UMAP", pcs=pcs)
    print(markers.keys())

def scvis_by_cluster_markers(rdata, tenx_analysis, prefix, pcs, embedding_file):
    try:
        tenx = TenxAnalysis(tenx_analysis)
        sce = SingleCellExperiment.fromRData(rdata)
        cluster_labels = tenx.markers_by_clusters(sce, rep="SCVIS", pcs=pcs, embedding_file=embedding_file)
    except Exception as e:
        return

def pca_by_cluster(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce, pcs=pcs)
    tsne_dims = sce.getReducedDims("PCA",n=pcs)
    barcodes = sce.colData["Barcode"]
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode, cluster in cluster_labels.items():
        try:
            x_val = x_coded[barcode]
            y_val = y_coded[barcode]
        except Exception as e:
            continue
        x.append(x_val)
        y.append(y_val)
        clusters.append("Cluster {}".format(cluster))
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("PCA - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/pca_by_cluster.png")


def scvis_by_cluster(rdata, tenx_analysis, prefix, pcs, embedding_file):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce, pcs=pcs)
    rows = open(embedding_file,"r").read().splitlines()
    dims = []
    rows.pop(0)
    for row in rows:
        row = row.split("\t")
        row = list(map(float,row[1:]))
        dims.append(row)
    barcodes = sce.colData["Barcode"]
    print(dims)
    x = []
    y = []
    clusters = []
    for barcode, dim in zip(barcodes,dims):
        x.append(dim[0])
        y.append(dim[1])
        clusters.append("Cluster {}".format(cluster_labels[barcode]))
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("SCVIS - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/svis_by_cluster_{}.png".format(prefix))


def umap_by_cluster(rdata, tenx_analysis, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce, pcs=pcs)
    tsne_dims = sce.reducedDims["UMAP"]
    barcodes = sce.colData["Barcode"]
    tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode, cluster in cluster_labels.items():
        clusters.append("Cluster {}".format(cluster))
        x.append(x_coded[barcode])
        y.append(y_coded[barcode])
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("PCA - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/umap_by_cluster.png")

def umap_by_gene(rdata, gene, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    tsne_dims = sce.reducedDims["UMAP"]
    barcodes = sce.colData["Barcode"]
    transcripts = sce.rowData["Symbol"]
    adata = tenx.create_scanpy_adata(barcodes=barcodes,transcripts=symbols)
    assert len(barcodes) == len(adata[:,gene])
    expression = dict(zip(barcodes,adata[:,gene]))
    tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode in barcodes:
        clusters.append(float(expression[barcode]))
        x.append(x_coded[barcode])
        y.append(y_coded[barcode])
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("PCA - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/umap_by_{}.png".format(gene))

def plot_by_genes(rdata, tenx_analysis, genes, prefix, rep, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    tsne_dims = sce.getReducedDims(rep)
    barcodes = sce.colData["Barcode"]
    transcripts = sce.rowData["Symbol"]
    adata = tenx.create_scanpy_adata(barcodes=barcodes,transcripts=transcripts)
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    if not os.path.exists("figures/expression"):
        os.makedirs("figures/expression")
    x = []
    y = []
    for barcode in barcodes:
        x.append(x_coded[barcode])
        y.append(y_coded[barcode])
    for gene in genes:
        expression = []
        for barcode in barcodes:
            val = adata[barcode,gene].X
            expression.append(float(val))
        f, ax = plt.subplots(figsize=(10,8))
        sns.scatterplot(x=x, y=y, hue=expression, alpha=0.85)
        ax.set_title("{} Counts".format(gene))
        ax.legend()
        plt.tight_layout()
        plt.savefig("figures/expression/expression_{}.png".format(gene))

def plot_by_markers(rdata, tenx_analysis, genes, prefix, rep, pcs, embedding_file=None, k = 12):
    genes = list(genes[:k])
    sce = SingleCellExperiment.fromRData(rdata)
    counts = sce.assays["logcounts"].toarray()
    tenx = TenxAnalysis(tenx_analysis)
    if rep == "SCVIS":
        tsne_dims = tenx.get_scvis_dimensions(embedding_file)
    else:
        tsne_dims = sce.getReducedDims(rep)
    all_genes = tenx.get_genes(sce)
    barcodes = sce.colData["Barcode"]
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    if not os.path.exists("figures/expression"):
        os.makedirs("figures/expression")
    x = []
    y = []
    f = plt.figure()
    drop_rows = []
    for i, barcode in enumerate(barcodes):
        try:
            x_val = x_coded[barcode]
            y_val = y_coded[barcode]
            x.append(x_val)
            y.append(y_val)
        except Exception as e:
            drop_rows.append(i)
            continue
        print("Barcode {}".format(barcode))
    scale = preprocessing.MinMaxScaler()
    counts = scale.fit_transform(counts)
    for i, gene in enumerate(genes):
        expression = counts[all_genes.index(gene)]
        plt.subplot(3, 4, i+1)
        expression = scale.fit_transform(numpy.array(expression).reshape(-1,1))
        _expression = list(expression.flatten())
        expression = []
        for i, row in enumerate(_expression):
            if i not in drop_rows:
                expression.append(row)
        print(len(expression))
        print(len(barcodes))
        print(len(x))
        print(len(y))
        g = sns.scatterplot(x=x, y=y, hue=expression, palette="RdYlBu_r", alpha=0.7, legend=False, s=4)
        g.set(xticklabels=[])
        g.set(yticklabels=[])
        plt.title(gene)
        print("Gene {}".format(gene))
    plt.tight_layout()
    plt.savefig("figures/expression/{}_counts_{}.png".format(prefix.replace(" ","_").lower(), rep.lower()))

def pca_by_gene(rdata, genes, prefix, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    tsne_dims = sce.getReducedDims("PCA",n=pcs)
    barcodes = sce.colData["Barcode"]
    transcripts = sce.rowData["Symbol"]
    adata = tenx.create_scanpy_adata(barcodes=barcodes,transcripts=symbols)
    assert len(barcodes) == len(adata[:,gene])
    expression = dict(zip(barcodes,adata[:,gene]))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    clusters = []
    for barcode in barcodes:
        clusters.append(float(expression[barcode]))
        x.append(x_coded[barcode])
        y.append(y_coded[barcode])
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85)
    ax.set_title("PCA - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/umap_by_{}.png".format(gene))



def cell_type_by_cluster(rdata, cell_assign_fit, tenx_analysis, prefix):
    tenx = TenxAnalysis(tenx_analysis)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cell_types = dict(zip(fit["Barcode"],fit["cell_type"]))
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce)
    clusters = dict(zip(sce.colData["Barcode"], cluster_labels))
    data_by_cluster = collections.defaultdict(list)
    data_by_celltype = collections.defaultdict(list)
    cluster = []
    cell_type = []
    for barcode, cell in cell_types.items():
        try:
            cluster.append(str(clusters[barcode]))
            cell_type.append(cell)
            data_by_celltype[cell] = str(clusters[barcode])
            data_by_cluster[str(clusters[barcode])] = cell
        except Exception as e:
            continue
    f, ax = plt.subplots(figsize=(16,8))
    counts = collections.defaultdict(lambda : collections.defaultdict(int))
    for cluster, ctype in zip(cluster, cell_type):
        counts[cluster][ctype] += 1
    fclusters = []
    fcelltypes = []
    fpercentages = []
    for cluster, ctype in counts.items():
        total = float(sum(ctype.values()))
        for cell in cell_type:
            fcelltypes.append(cell)
            fclusters.append(cluster)
            if cell in ctype:
                fpercentages.append(float(ctype[cell]) / total)
            else:
                fpercentages.append(0.0)
    df = pandas.DataFrame({"Cluster":fclusters, "Cell Type": fcelltypes, "Percentage": fpercentages})
    ax = sns.barplot(x="Cluster", y= "Percentage", hue="Cell Type", data=df, palette="tab10")
    ax.set_title("Cell Type by Cluster - {}".format(prefix))
    plt.tight_layout()
    plt.savefig("figures/cell_type_by_cluster.png")


def marker_analysis(sce, tenx, rho, cell_assign_fit, figure):
    sce = SingleCellExperiment.fromRData(sce)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    gene_markers = []
    for markers in rho.values():
        print(markers)
        gene_markers += markers
    _marker_genes = list(set(gene_markers))
    convert = tenx.gene_map(sce)
    marker_genes = []
    for gene in _marker_genes:
        try:
            marker_genes.append(convert[gene])
        except KeyError:
            marker_genes.append(gene)
            print('No conversion for ',gene)
    print(marker_genes)
    adata = tenx.create_scanpy_adata(sce, fast_load=False)
    print(len(adata.obs.index))
    print(len(fit["cell_type"]))
    cell_types = []
    _cell_types = dict(zip(fit["Barcode"],fit["cell_type"]))
    for barcode in adata.obs.index:
        try:
            cell_types.append(_cell_types[barcode])
        except KeyError as e:
            cell_types.append("Other")
    adata.obs["Cell Type"] = cell_types
    print(len(cell_types))
    marker_genes = list(set(marker_genes).intersection(set(adata.var.index)))
    print(len(marker_genes))
    print(marker_genes)
    sc.pl.dotplot(adata, marker_genes, groupby='Cell Type', save="matrix.png")
    sc.pl.stacked_violin(adata, marker_genes, groupby='Cell Type', rotation=90, save="vin_stacked.png")
    return ["dot_plot.png", "vin_stacked.png"]


if __name__ == '__main__':
    import glob
    pngs = glob.glob("/home/nceglia/jobs/SC_723/*/*0.png")
    combine_figures(pngs,"SC_723_vis.pdf")
