import argparse
import yaml
import sys
import os

from utils.config import write_config

parser = argparse.ArgumentParser()

parser.add_argument('--sampleid', type=str, help='Sample ID linked to fastqs in scrnadata.')
parser.add_argument("--build", type=str, help="Reference build.", default="GRCh38")
parser.add_argument("--jobpath", type=str, help="Path to job folder.", default = "/results")
parser.add_argument("--datapath", type=str, help="Path to directory holding fastqs.", default ="/data")
parser.add_argument("--referencepath", type=str, help="Path to reference directory with reference in dir labelled by build.", default="/reference")
parser.add_argument("--cellranger", type=str, help="Path to cellranger bin.", default="/codebase/cellranger-3.0.2/cellranger-cs/3.0.2/bin/")
parser.add_argument("--lsf", action='store_true', help="Path to cellranger bin.")
parser.add_argument("--markers", type=str, help="Path to marker matrix.")

args = parser.parse_args()

if not os.path.exists(os.path.join(args.datapath, args.sampleid)):
    os.makedirs(os.path.join(args.datapath, args.sampleid))


write_config(args.sampleid, args.build, args.jobpath, os.path.join(args.datapath, args.sampleid), args.referencepath, args.cellranger, args.lsf, args.markers)
print("Config written to settings.yaml")
