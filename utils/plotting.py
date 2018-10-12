import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import math

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
