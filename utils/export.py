import os

from interface.singlecellexperiment import SingleCellExperiment

def embed_cellranger(directory):
    summary = os.path.join(directory,"web_summary.html")
    assert os.path.exists(summary)
    code = """
    ```
    \n\n\n
    {{r, echo=FALSE}}
    htmltools::tags$iframe(title = "Cellranger Output", src = "{web}")
    ```
    \n\n\n
    """.format(web=summary)

def codeblock(code, handle):
    handle.write("\n\n\n```\n")
    handle.write("library(scater)")
    handle.write("")
    handle.write(code)
    handle.write("\n```\n\n\n")

def exportRMD(directory, prefix, fastq, sce_qc):
    rmd = os.path.join(directory,"{}_analysis.rmd".format(prefix))
    output = open(rmd,"w")
    output.write("# {} Analysis".format(prefix))
    output.write("\n\nPipeline Specifics\n")
    output.write(" - FastQ location: {}".format(fastq.results))
    embed_cellranger(fastq.results)

    # code = "plotHighestExprs({}, exprs_values = 'counts')".format(sce_qc)
    # codeblock(code, output)

    # code = "plotExprsFreqVsMean({})".format(sce_qc)
    # codeblock(code, output)

    output.close()



if __name__ == '__main__':
    pass
