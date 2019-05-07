import json
import sys
data = json.loads(open("patient_data.json","r").read())
fields = data.keys()
samples = list(filter(lambda x: "IGO" in x, fields))
right_ovary = data["LA-CD45N_IGO_09443_D_10"]
right_ovary_cells = right_ovary["celldata"]["Barcode"] # All barcodes for sample
right_ovary_genes = right_ovary["genedata"]["Symbol"] # All genes for sample
barcode = right_ovary_cells[10]
gene = right_ovary_genes[100]
print(barcode,gene)

stats_by_samples = data["statistics"]
right_ovary_stats = stats_by_samples["LA-CD45N_IGO_09443_D_10"]
cellranger_kit = right_ovary_stats["Chemistry"]
estimated_number_of_cells = int(right_ovary_stats["Estimated Number of Cells"])
print("Cells",estimated_number_of_cells)
cells_retained_at_default_value_for_qc = int(right_ovary_stats["Mito10"])
assert estimated_number_of_cells > cells_retained_at_default_value_for_qc, "We keep less cells than cellranger finds."
saturation = right_ovary_stats["Sequencing Saturation"]
print("Saturation",saturation)

rkeys = []
new_rho = dict()
for ct, markers in data["rho"].items():
    celltype = " ".join([x.capitalize() for x in ct.split()])
    new_rho[celltype] = markers

data["rho"] = new_rho

rho = data["rho"]
cell_types = list(rho.keys())
all_marker_genes = list(rho.values())
print(cell_types)
sys.stdout.flush()
cell_type = cell_types[0]
plasma_cell_marker_genes = rho[cell_type]

matrix = right_ovary["matrix"]
gene = list(matrix[barcode].keys())[100]
count = matrix[barcode][gene]
try:
  count = matrix[barcode]["MIR1302-10"] #No expression
except KeyError:
  count = 0


for sample in samples:
    for bc, cell_type in data[sample]["cellassign"].items():
        ct = " ".join([x.capitalize() for x in cell_type.split(".")])
        data[sample]["cellassign"][bc] = ct


tsne = right_ovary["tsne"]
x_coordinate, y_coordinate = tsne[barcode]

celltype = right_ovary["cellassign"][barcode]

assert celltype in cell_types, "This celltype should be expected in the rho matrix."
# Get the marker gene expression in this cell for the given cell type.
marker_genes = rho[celltype]
counts_by_gene = matrix[barcode]
expression_subset = dict((gene, counts_by_gene[gene]) for gene in marker_genes if gene in counts_by_gene)


cluster = right_ovary["clusters"][barcode]

celldata = right_ovary["celldata"]
total_counts_to_be_binned = celldata["total_counts"] #make me a histo or a violin and compare across sites :)
total_counts_to_be_binned = celldata["total_features_by_counts"]
total_counts_mito_to_be_binned = celldata["total_counts_mito"]
pct_counts_mito_to_be_binned = celldata["pct_counts_mito"] #Use me to show finer resolution on whether 10% is a good QC :)
total_counts_ribo_to_be_binned = celldata["total_counts_ribo"]
total_counts_ribo_to_be_binned = celldata["total_counts_ribo"]

patient_data_str = json.dumps(data)
output = open("patient_data2.json","w")
output.write(str(patient_data_str))
output.close()
