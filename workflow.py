import pypeliner.workflow
import argparse
from utils import converter

def create_workflow():
    workflow = pypeliner.workflow.Workflow()
    return workflow


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)
    args = vars(argparser.parse_args())
    workflow = create_workflow()
    pyp = pypeliner.app.Pypeline(config=args)
    pyp.run(workflow)
