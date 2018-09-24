import subprocess
import os
import random
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis

class CellRanger(object):

    def __init__(self):
        pass

    @staticmethod
    def cmd(command, args):
        cmd = ["cellranger", command]
        for flag, value in args.items():
            cmd.append("--{}".format(flag))
            cmd.append(value)
        return cmd

    @staticmethod
    def mkfastq(bcl_object):
        args = dict()
        args["id"] = bcl_object.id
        args["run"] = bcl_object.path
        args["csv"] = bcl_object.csv
        cmd = CellRanger.cmd("mkfastq",args)
        subprocess.call(cmd)
        return FastQDirectory("/home/ceglian/data/fastqs")


    @staticmethod
    def count(fastq_object):
        args = dict()
        args["id"] = fastq_object.id
        args["fastqs"] = fastq_object.path
        args["sample"] = "hgmm_100"
        args["transcriptome"] = "/juno/work/shah/ceglian/refdata/hg19_and_mm10"
        args["lanes"] = "1"
        args["chemistry"] = "SC3P_auto"
        cmd = CellRanger.cmd("count",args)
        subprocess.call(cmd)
        return TenxAnalysis(fastq_object.out())



    def aggr(self):
        pass
