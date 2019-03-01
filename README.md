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
