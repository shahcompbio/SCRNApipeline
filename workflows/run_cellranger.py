import pypeliner.workflow
import pypeliner.app
import pypeliner.managed

import sys
import os

from software.cellranger import CellRanger

def RunCellranger(fastqs, workflow):
    assert len(fastqs) > 0, "Need to provide atleast 1 fastq directory."
    completed = os.path.join(fastqs[0].output, "cellranger.complete")
    workflow.transform (
        name = "{}_cellranger_counts".format(fastqs[0].id),
        func = CellRanger.count,
        args = (
            fastqs,
        )
    )
    return workflow
