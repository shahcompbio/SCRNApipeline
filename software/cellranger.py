import subprocess
import os
import random
from interface.fastqdirectory import FastQDirectory
from interface.tenxanalysis import TenxAnalysis

from utils.config import Configuration

config = Configuration()

class CellRanger(object):

    @staticmethod
    def cmd(command, args):
        cmd = [os.path.join(config.cellranger,"cellranger"), command]
        for flag, value in args.items():
            cmd.append("--{}".format(flag)
            )
            cmd.append(value)
        cmd.append("--disable-ui")
        cmd.append("--jobmode=lsf")
        return cmd

    @staticmethod
    def mro(fastqs, mrofile):
        sample_json = dict()
        content = """
        @include "sc_rna_counter.cs.mro"
        call SC_RNA_COUNTER_CS(
            sample_id = {sample_id},
            sample_def = [
                {sample_definitions}
            ],
            sample_desc = "",
            reference_path = {reference},
            recovered_cells = null,
            force_cells = null,
            no_secondary_analysis = false,
        )"""
        sample_definitions = list()
        for fastq in fastqs:
            definition = dict()
            defintion["fastq_mode"] = "ILMN_BCL2FASTQ"
            definition["gem_group"] = "null",
            definition["lanes"] = "null",
            definition["read_path"] = fastq.path
            definition["sample_indices"] = '["any"]'
            definition["sample_names"] = '[" ' + fastq.sampleid + ' "]'
            sample_definitions.append(definition)
        output = open(mrofile, "w")
        output.write(content.format(sample_id=run_id, sample_definitions=sample_definitions, reference=config.reference))
        output.close()


    @staticmethod
    def mkfastq(bcl_object):
        args = dict()
        args["id"] = bcl_object.id
        args["run"] = bcl_object.path
        args["csv"] = bcl_object.csv
        cmd = CellRanger.cmd("mkfastq",args)
        subprocess.call(cmd)
        return FastQDirectory(bcl_object.out())

    @staticmethod
    def aggr(csv, prefix):
        args = dict()
        args["id"] = "run_{}".format(prefix)
        args["csv"] = csv
        args["normalize"] = "none"
        cmd = CellRanger.cmd("aggr",args)
        subprocess.call(cmd)

    @staticmethod
    def count(fastqs):
        args = dict()
        fastq_files = []
        args["id"] = "_".join([fastq.id for fastq in fastqs])
        paths = [fastq.path for fastq in fastqs]
        args["fastqs"] = ",".join(paths)
        try:
            args["sample"] = ",".join([fastq.samples.sampleid[0] for fastq in fastqs])
        except Exception as e:
            pass
        args["transcriptome"] = config.reference
        print(config.reference)
        if config.chemistry is not None:
            args["chemistry"] = config.chemistry
        cmd = CellRanger.cmd("count",args)
        print("Running ", " ".join(cmd))
        subprocess.call(cmd)

    @staticmethod
    def reanalyze(tenx_object):
        args = dict()
        args["id"] = tenx_object.id
        args["matrix"] = tenx_object.matrix
        args["params"] = tenx_object.params
        cmd = CellRanger.cmd("reanalyze", args)
        subprocess.call(cmd)
