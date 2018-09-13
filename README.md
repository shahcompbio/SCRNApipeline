# SCRNApipeline #
Single Cell RNA-seq Pipeline

**Requirements**
1. Python 3.x
2. rpy2
3. pandas
4. development branch of Pypeliner

## Modules ##
### Utils ###
1. export.py -> rdata to fl

### Interface ###
1. dropletutils
2. singlecellexperiment

### Software ###
1. cellassign
2. clonealign


## Examples ##
Loading RData
```
sce = SingleCellExperiment.fromRData("~/data/example_sce.RData")

```

Accessing SingleCellExperiment Assays as Pandas Dataframe (with row and column indexing)
```
for assay in sce.assayNames:
    dataframe = sce.assay(assay)
```

Accessing SingleCellExperiment Fields
```
print(sce.rowData)
print(sce.colData)
print(sce.assays)
print(sce.reducedDims)
print(sce.sizeFactors)
```


Exporting RData Example
```
input_rdata = "/home/ceglian/data/example_sce.RData"
output_directory = "/home/ceglian/data/example_sce_data/"
exportRData(input_rdata, output_directory)
```


Running DropletUtils and Loading RS4
```
tenx = DropletUtils()
res = tenx.read10xCounts("~/data/raw_gene_bc_matrices/GRCh38/")
sce = SingleCellExperiment.fromRS4(res)
```
