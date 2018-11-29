from interface.singlecellexperiment import SingleCellExperiment
from interface.tenxanalysis import TenxAnalysis
import scanpy.api as sc




rows = open("/home/nceglia/codebase/refdata/marker_genes.tsv","r").read().splitlines()
header = rows.pop(0).split("\t")
import collections
matrix = collections.defaultdict(list)
for row in rows:
    row = dict(zip(header,row.split("\t")))
    if "Hs" in row["species"]:
        if float(row["ubiquitousness index"]):
            matrix[row["cell type"]].append(row["official gene symbol"])

print(dict(matrix))

#
# for cell_type in matrix.keys():
#     print(cell_type, matrix[cell_type])
# print(len(matrix.keys()))
#
#
# analysis = TenxAnalysis("/home/nceglia/jobs/SI_GA_A/run_SI_GA_A/outs")
# adata = analysis.create_scanpy_adata()
# adata.var_names_make_unique()
# sc.pp.recipe_zheng17(adata)
# print(adata.var_names)
# print(len(adata.obs_names))
# expressed_genes = collections.defaultdict(list)
# for cell_type in matrix.keys():
#     count = 0
#     for gene in matrix[cell_type]:
#         if gene in adata.var_names:
#             count += 1
#             expressed_genes[gene].append(cell_type)
#     if count > 0:
#
#         print(cell_type, count)
#
# for gene, cell_types in expressed_genes.items():
#     print(gene, cell_types)


# sce = SingleCellExperiment.fromRData("/home/nceglia/jobs/SI_GA_A/sce_final.rdata")
# print(len(list(sce.rowData["Symbol"])
