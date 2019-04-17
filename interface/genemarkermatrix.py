import collections
import operator
import pandas
import numpy
import os

class GeneMarkerMatrix(object):

    def __init__(self, marker_list):
        self.marker_list = marker_list
        self.genes = list(set([gene for genelist in self.marker_list.values() for gene in genelist]))
        self.cells = list(self.marker_list.keys())

    @classmethod
    def read_yaml(cls, yaml):
        marker_list = collections.defaultdict(list)
        yaml_rows = open(yaml,"r").read().splitlines()
        celltype = yaml_rows.pop(0).replace(":","").strip()
        for line in yaml_rows:
            if line.strip().startswith("#"):
                continue
            elif line.strip().startswith("-"):
                line = line.split("#")[0]
                marker_list[celltype].append(line.replace("-","").strip())
            elif ":" in line:
                celltype = line.replace(":","").strip()
        return cls(marker_list)

    def write_matrix(self, filename, include_other=True):
        self.indicator = []
        self.genes = list(set([gene for genelist in self.marker_list.values() for gene in genelist]))
        self.cells = list(self.marker_list.keys())
        if include_other:
            self.cells.append("Other")
        for cell in self.cells:
            marker_genes = []
            for gene in self.genes:
                if cell == "Other":
                    marker_genes.append("0")
                elif gene in self.marker_list[cell]:
                    marker_genes.append("1")
                else:
                    marker_genes.append("0")
            self.indicator.append(marker_genes)
        self.indicator = zip(*self.indicator)
        output = open(filename, "w")
        output.write(",".join([""] + self.cells) + "\n")
        for marker, binary_indicators in zip(self.genes, self.indicator):
            row = [marker] + list(binary_indicators)
            output.write(",".join(row)+"\n")
        output.close()

    def celltypes(self):
        return self.cells
