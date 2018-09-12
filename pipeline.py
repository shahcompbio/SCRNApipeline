import pypeliner.workflow
import rpy2.robjects as robjects
import argparse
from utils import converter


def run_cellalign():
    r_source = robjects.r['source']
    r_source("")

def create_workflow():

    workflow = pypeliner.workflow.Workflow()
    #
    # workflow.commandline(
    #
    # ),
    #
    # workflow.commandline(
    #     name="cellassign",
    # ),



    # workflow.transform(
    #     name='split_fastq',
    #     func=split_file_byline,
    #     args=(
    #         pypeliner.managed.InputFile('/nicky/input.fastq'),
    #         100,
    #         pypeliner.managed.TempOutputFile('input.fastq', 'reads'),
    #     )
    # )
    #
    # workflow.commandline(
    #     name='bwa_aln',
    #     axes=('reads',),
    #     args=(
    #         'bwa', 'aln',
    #         pypeliner.managed.InputFile("/nicky/genome.fa"),
    #         pypeliner.managed.TempInputFile('input.fastq', 'reads'),
    #         '>',
    #         pypeliner.managed.TempOutputFile('tmp.sai','reads'),
    #     )
    # )
    #
    #
    # workflow.commandline(
    #     name='bwa_samse',
    #     axes=('reads',),
    #     args=(
    #         'bwa', 'samse',
    #         pypeliner.managed.InputFile('/nicky/genome.fa'),
    #         pypeliner.managed.TempInputFile('tmp.sai','reads'),
    #         pypeliner.managed.TempInputFile('input.fastq','reads'),
    #         '>',
    #         pypeliner.managed.TempOutputFile('output.sam','reads'),
    #     )
    # )
    #
    # workflow.transform(
    #     name='merge_sam',
    #     func=merge_file_byline,
    #     args=(
    #         pypeliner.managed.TempInputFile('output.sam', 'reads'),
    #         pypeliner.managed.OutputFile('output.sam'),
    #     )
    # )
    #
    # workflow.transform(
    #     name='filter_unmapped',
    #     func=filter_unmapped,
    #     args=(
    #         pypeliner.managed.InputFile('output.sam'),
    #         pypeliner.managed.OutputFile('aligned.sam'),
    #     )
    # )
    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)
    args = vars(argparser.parse_args())
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
