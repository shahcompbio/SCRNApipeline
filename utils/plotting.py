import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import math
from utils.config import *
import pickle
from interface.genemarkermatrix import GeneMarkerMatrix
from interface.singlecellexperiment import SingleCellExperiment
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import collections
sns.set(style="darkgrid")

def scatter(x, y, xlabel, ylabel, hlines, filename):
    df = dict()
    df[xlabel] = x
    df[ylabel] = y
    df = pandas.DataFrame.from_dict(df)
    f, ax = plt.subplots(figsize=(8,10))
    ax.set(xscale="log", yscale="log")
    sns.scatterplot(x=xlabel, y=ylabel, data = df)
    colors = ["r","g","b","c"]
    for label, hline in hlines.items():
        ax.axhline(y=hline[0],label=label, color=colors.pop(0))
    ax.set_title("Barcode Ranks")
    ax.legend()
    plt.savefig(filename)

def celltypes(data, filename, cell_labels):
    f, ax = plt.subplots(figsize=(10,8))
    sns.countplot(data)
    ax.set_title("Cell Type Assignments")
    sns.countplot(data)
    ax.set_xticklabels(labels=cell_labels,rotation=30)
    plt.tight_layout()
    plt.savefig(filename)

def tsne_scanpy(tenx, cell_assign_fit, perplexity, filename):
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
        print("Coordinate")
        try:
            labels.append(assignments[barcode])
            x.append(float(coordinate[0]))
            y.append(float(coordinate[1]))
        except Exception as e:
            #labels.append("Other")
            continue
    f, ax = plt.subplots(figsize=(10,8))
    sns.scatterplot(x=x, y=y, hue=labels,alpha=0.3)
    ax.set_title("TSNE - Cell Types")
    ax.legend()
    plt.tight_layout()
    plt.savefig(filename)

def tsne_scran(rdata, cell_assign_fit):
    sce = SingleCellExperiment.fromRData(rdata)
    print(type(sce.reducedDims))
    dims = pandas2ri.ri2py(sce.reducedDims)
    dims = dict(zip(dims.names, list(dims)))
    for key, val in dims.items():
        print(key)
        print(type(val))
        print(key, val)

def umap(rdata, tenx, cell_assign_fit, min_dist, filename):
    tenx.create_scanpy_adata()
    projection = tenx.umap(min_dist=min_dist)
    fit = pickle.load(open(cell_assign_fit,"rb"))
    cells = fit["cell_type"]
    barcodes = fit["Barcode"]
    assert len(cells) == len(barcodes), "No matching barcodes {} <-> {}".format(len(cells), len(barcodes))
    assignments = dict(zip(barcodes, cells))
    x = []
    y = []
    labels = []
    for barcode, coordinate in projection.items():
        print("Coordinate")
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
