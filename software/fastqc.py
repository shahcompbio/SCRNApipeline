from interface.qcreport import QCReport
import subprocess
import glob
import os

class FastQC(object):

    @staticmethod
    def cmd(fastqs, args):
        for fastq in fastqs:
            cmd = ["fastqc"]
            cmd += fastqs
            for flag, value in args.items():
                cmd.append("--{}".format(flag))
                cmd.append(value)
            yield cmd

    @staticmethod
    def run(fastq_object, output):
        args = dict()
        args["outdir"] = output
        fastqs = []
        fastqs = glob.glob(os.path.join(fastq_object.path,"*.fastq*"))
        cmds = FastQC.cmd(fastqs,args)
        for cmd in cmds:
            subprocess.call(cmd)
        return QCReport(output)
