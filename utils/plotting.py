import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import math
from utils.config import *
import pickle
from interface.genemarkermatrix import GeneMarkerMatrix
from interface.singlecellexperiment import SingleCellExperiment
from interface.tenxanalysis import TenxAnalysis
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
    ax.set_xticklabels(labels=list(sorted(order)),rotation=30)
    plt.tight_layout()
    figure = os.path.join(output, "cell_types.png")
    plt.savefig(figure)

def tsne_perplex(tenx, cell_assign_fit, perplexity, filename):
    tenx.create_scanpy_adata()
    print("Plotting TSNE")
    projection = tenx.tsne(perplexity=perplexity)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cells = fit["cell_type"]
    barcodes = fit["Barcode"]
    assert len(cells) == len(barcodes), "No matching barcodes {} <-> {}".format(len(cells), len(barcodes))
    assignments = dict(zip(barcodes, cells))
    x = []
    y = []
    labels = []
    for barcode, coordinate in projection.items():
        try:
            labels.append(assignments[barcode])
            x.append(float(coordinate[0]))
            y.append(float(coordinate[1]))
        except Exception as e:
            continue
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=labels,alpha=0.3)
    ax.set_title("TSNE - Cell Types")
    ax.legend()
    plt.tight_layout()
    plt.savefig(filename)

def tsne_by_cell_type(rdata, cell_assign_fit, prefix):
    sce = SingleCellExperiment.fromRData(rdata)
    tsne_dims = sce.reducedDims["TSNE"]
    barcodes = sce.colData["Barcode"]
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cell_types = dict(zip(fit["Barcode"][:len(barcodes)],fit["cell_type"][:len(barcodes)]))
    tsne_dims = numpy.array(tsne_dims).reshape(2, len(barcodes))
    x_coded = dict(zip(barcodes, tsne_dims[0]))
    y_coded = dict(zip(barcodes, tsne_dims[1]))
    x = []
    y = []
    for barcode, cell_type in cell_types.items():
        try:
            x_val = x_coded[barcode]
            y_val = y_coded[barcode]
            x.append(x_val)
            y.append(y_val)
        except Exception as e:
            continue
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=fit["cell_type"][:len(x)],alpha=0.85, palette="tab10")
    ax.set_title("TSNE - Cell Types - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("tsne_by_cell_type.png")

def tsne_by_cluster(rdata, tenx_analysis, prefix):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce)
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
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85, palette="viridis")
    ax.set_title("TSNE - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("tsne_by_cluster_{}.png".format(prefix))

def scvis_by_cluster(rdata, tenx_analysis, prefix, embedding_file, pcs):
    tenx = TenxAnalysis(tenx_analysis)
    sce = SingleCellExperiment.fromRData(rdata)
    cluster_labels = tenx.clusters(sce,embedding_file=embedding_file, pcs=pcs)
    tsne_dims = []
    print("FILE", embedding_file)
    rows = open(embedding_file,"r").read().splitlines()
    rows.pop(0)
    for row in rows:
        row = row.split("\t")
        row = list(map(float,row[1:]))
        tsne_dims.append(row)
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
    sns.scatterplot(x=x, y=y, hue=clusters, alpha=0.85, palette="viridis")
    ax.set_title("SCVIS - Clusters - {}".format(prefix))
    ax.legend()
    plt.tight_layout()
    plt.savefig("svis_by_cluster_{}.png".format(prefix))

def cell_type_by_cluster(rdata, cell_assign_fit, prefix):
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cell_types = dict(zip(fit["Barcode"],fit["cell_type"]))
    sce = SingleCellExperiment.fromRData(rdata)
    clusters = dict(zip(sce.colData["Barcode"], sce.colData["Cluster"]))
    data_by_cluster = collections.defaultdict(list)
    data_by_celltype = collections.defaultdict(list)
    cluster = []
    cell_type = []
    for barcode, cell in cell_types.items():
        cluster.append(str(clusters[barcode]))
        cell_type.append(cell)
        data_by_celltype[cell] = str(clusters[barcode])
        data_by_cluster[str(clusters[barcode])] = cell
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
    plt.savefig("cell_type_by_cluster.png")

def umap(rdata, tenx, cell_assign_fit, filename):
    tenx.create_scanpy_adata()
    projection = tenx.umap()
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cells = fit["cell_type"]
    barcodes = fit["Barcode"]
    assert len(cells) == len(barcodes), "No matching barcodes {} <-> {}".format(len(cells), len(barcodes))
    assignments = dict(zip(barcodes, cells))
    x = []
    y = []
    labels = []
    for barcode, coordinate in projection.items():
        try:
            labels.append(assignments[barcode])
            x.append(float(coordinate[0]))
            y.append(float(coordinate[1]))
        except Exception as e:
            continue
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=labels,alpha=0.3)
    ax.set_title("UMAP - Cell Types")
    ax.legend()
    plt.tight_layout()
    plt.savefig(filename)


if __name__ == '__main__':
    import glob
    pngs = glob.glob("/home/nceglia/jobs/SC_723/*/*0.png")
    combine_figures(pngs,"SC_723_vis.pdf")
