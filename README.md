# Single Cell RNA-Seq Pipeline #

Pipeline for running single cell rna-seq experiments.
This is primarily built on Cell Ranger with additionaly analysis from CellAssign, CloneAlign, and SCViz tools.
The workflow is inferred based on the inclusion (or omission) of command line arguments.

# Running with Docker #

Create a job directory with the directory structure including jobs/ data/ reference/.
```
mkdir somedirectory/test-run
mkdir somedirectory/test-run/jobs
mkdir somedirectory/test-run/data
mkdir somedirectory/test-run/reference
```

Copy dockerfile and settings yaml into job directory.
```
cp docker/settings.yaml somedirectory/test-run/
cp docker/dockerfile somedirectory/test-run/
```

Edit settings.yaml with prefix, genome build, and datapath.
```
prefix: "test"
build: "test"
datapath: "/data/test"
```

Build image.
```
docker build .
```

Run docker.
```
docker run -v /Users/ceglian/Development/docker/data:/data -v /Users/ceglian/Development/docker/jobs/:jobs -v /Users/ceglian/Development/docker/reference:/reference
```


Main results are stored in mounted volumes.
```
ls PREFIX_report/PREFIX_results.html
```




Notes:
If the build specified in the settings.yaml is not available, it will be pulled and built.  This will take a long time.



## Examples ##

Running pipeline starting from binary base call directory (BCL) to report generation.
The top level BCL directory must include a single csv worksheet with minimum columns for Lane, Sample, Index (10x index used for library construction).

```
    $ python3 pipeline.py --bcl tests/cellranger-tiny-bcl-1.2.0/
```

The FastsQ directory can be the output of `cellranger mkfastq` or a directory where fastq files are named '{sample}_S{sample_num}_L00{lane}_{R{read} || I1}_001'.
If this argument is omitted, but the `bcl` argument is included, the fastq path will be inferred.

```
    $ python3 pipeline.py --fastq tests/tiny-fastqs-mk/outs/fastq_path/
```

The tenx analysis folder can be the output of `cellranger count` or a directory that includes a {filtered || raw}_gene_bc_matrices folder.
If this argument is omitted, but bcl or fastq is included, the path will be inferred.
If this directory includes Cell Ranger analysis, this will be included in the report generation.

```
    $ python3 pipeline.py --tenx tests/tiny-fastqs-count/outs/
```

Single Cell Experiment objects can be created from tenx analysis folders or loaded from serialized RData objects.
You can load these and run the pipeline downstream analysis starting from these serialized objects.

```
    $ python3 pipeline.py --rdata tests/example_sce.RData
```


# SCRNA Viz json

## Overview
Found in **scrnaviz** container in **scrnadata** storage account.
Right now, just an example. I will figure out how to index patient ID for subsequent runs.

- patient_data.json

```
import json
data = json.loads(open("patient_data.json","r").read())
```

- htmls
Just a directory of cellranger web summary htmls.

### Fields
```
fields = data.keys()
```

1. Sample Names (Sites)
Samples should start with Sample right now.
```
samples = list(filter(lambda x: "Sample" in x, fields))
right_ovary = data["Sample_ROVCD45POS_IGO_09443_B_1"]
right_ovary_cells = right_ovary["celldata"]["Barcode"] # All barcodes for sample
right_ovary_genes = right_ovary["genedata"]["Symbol"] # All genes for sample
```

2. Statistics
These values were taken from the web summary or qc R scripts.  They are currently strings :(.
They are all fairly important, try to compare all values across sites.
```
stats_by_samples = data["statistics"]
right_ovary_stats = stats_by_samples["Sample_ROVCD45POS_IGO_09443_B_1"]
cellranger_kit = right_ovary_stats["Chemistry"]
estimated_number_of_cells = int(right_ovary_stats["Estimated Number of Cells"])
cells_retained_at_default_value_for_qc = int(right_ovary_stats["Mito10"])
assert estimated_number_of_cells > cells_retained_at_default_value_for_qc, "We keep less cells than cellranger finds."
saturation = right_ovary_stats["Sequencing Saturation"]
```

3. Cellassign Rho
This includes all the expected cell types and the marker genes associated.
Will probably want to combine celltype with marker gene expression (heat).
```
rho = data["rho"]
cell_types = list(rho.keys())
all_marker_genes = list(rho.values())
plasma_cell_marker_genes = rho["Cytotoxic T cells"]
```

### Sample (Site) specific data

1. Matrix
Count expression (heat) matrix in cells by genes.
* Matrix is in a sparse format, if the cell does not have expression for a specific gene, will raise KeyError *
* Solve by putting into defaultdict(int) *

```
matrix = right_ovary["matrix"]
count = matrix["AAACCCAAGGCATGCA-1"]["CDK11A"]
try:
  count = matrix["AAACCCAAGGCATGCA-1"]["MIR1302-10"] #No expression
except KeyError:
  count = 0
```

2. tSNE
Coordinates for QC'd cells (Mito10).
```
tsne = right_ovary["tsne"]
x_coordinate, y_coordinate = tsne["AAACCCAAGGCATGCA-1"]
```

3. Cell Types
Cell types for QC'd cells (Mito10).
```
celltype = right_ovary["cellassign"]["AAACCCAAGGCATGCA-1"]
assert celltype in celltypes, "This celltype should be expected in the rho matrix."
# Get the marker gene expression in this cell for the given cell type.
marker_genes = rho[celltype]
counts_by_gene = matrix["AAACCCAAGGCATGCA-1"]
expression_subset = dict((gene, counts_by_gene[gene]) for gene in marker_genes if gene in counts_by_gene)
```

4. Clusters
Cluster assignments for QC'd cells from leiden algorithm.
```
cluster = right_ovary["clusters"]["AAACCCAAGGCATGCA-1"]
```

5. Cell Data

```
celldata = right_ovary["celldata"]
total_counts_to_be_binned = celldata["total_counts"] #make me a histo or a violin and compare across sites :)
total_counts_to_be_binned = celldata["total_features_by_counts"]
total_counts_mito_to_be_binned = celldata["total_counts_mito"]
pct_counts_mito_to_be_binned = celldata["pct_counts_mito"] #Use me to show finer resolution on whether 10% is a good QC :)
total_counts_ribo_to_be_binned = celldata["total_counts_ribo"]
```


Things to add:
1) Umap and scvis
2) ... will think more ...

Things to fix:
1) Create issues as needed to fix the format of the data.
