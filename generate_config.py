import argparse
import yaml
import sys

from utils.config import write_config

parser = argparse.ArgumentParser()

parser.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
parser.add_argument("--build", type=str, help="Reference build.", default="GRCh38")
parser.add_argument("--jobpath", type=str, help="Path to job folder.", default = "/results")
parser.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
parser.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
parser.add_argument("--cellranger", type=str, help="Path to cellranger bin.", default="/codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/")

args = parser.parse_args()

write_config(args.sampleid, args.build, args.jobpath, args.datapath, args.referencepath, args.cellranger)
print("Config written to settings.yaml")
