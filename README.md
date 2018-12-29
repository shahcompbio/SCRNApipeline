# Single Cell RNA-Seq Pipeline #

Pipeline for running single cell RNA-seq experiments.
This is primarily built on Cell Ranger with additional analysis from CellAssign, CloneAlign, and SCViz tools.
The workflow is inferred based on the inclusion (or omission) of command line arguments.

## Examples ##

Running pipeline starting from binary base call directory (BCL) to report generation.
The top level BCL directory must include a single CSV worksheet with minimum columns for Lane, Sample, Index (10x index used for library construction).

```
    $ python3 pipeline.py --bcl tests/cellranger-tiny-bcl-1.2.0/
```

The FASTQ directory can be the output of `cellranger mkfastq` or a directory where FASTQ files are named `{sample}_S{sample_num}_L00{lane}_{R{read} || I1}_001`.
If this argument is omitted, but the `bcl` argument is included, the FASTQ path will be inferred.

```
    $ python3 pipeline.py --fastq tests/tiny-fastqs-mk/outs/fastq_path/
```

The tenx analysis folder can be the output of `cellranger count` or a directory that includes a `{filtered || raw}_gene_bc_matrices` folder.
If this argument is omitted, but `bcl` or `fastq` is included, the path will be inferred.
If this directory includes Cell Ranger analysis, this will be included in the report generation.

```
    $ python3 pipeline.py --tenx tests/tiny-fastqs-count/outs/
```

Single Cell Experiment objects can be created from tenx analysis folders or loaded from serialized RData objects.
You can load these and run the pipeline downstream analysis starting from these serialized objects.

```
    $ python3 pipeline.py --rdata tests/example_sce.RData
```
