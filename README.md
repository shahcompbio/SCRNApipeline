# SCRNApipeline
Single Cell RNA-seq Pipeline


# # RData to TSV
Find here: **/utils/converter.py**
**Requirements**
1. Python 3.x
2. rpy2
3. tzlocal


```
python3 converter.py --rdata example_sce.RData --output example.tsv
```
or script it:
```
sce = SingleCellExperiment("example_sce.RData")
sce.dump("example.tsv")
```
Retrieve pandas dataframe:
```
sce = SingleCellExperiment("example_sce.RData")
dataframe = sce.matrix()
print(dataframe.head)
```
