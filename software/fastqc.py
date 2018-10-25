from interface.qcreport import QCReport
import subprocess
import glob
import os

class FastQC(object):

    @staticmethod
    def cmd(fastqs, args):
        for fastq in fastqs:
            cmd = ["fastqc"]
            cmd.append(fastq)
            for flag, value in args.items():
                cmd.append("--{}".format(flag))
                cmd.append(value)
            cmd.append("--noextract")
            yield cmd

    @staticmethod
    def run(fastq_object):
        args = dict()
        args["outdir"] = os.path.join(fastq_object.output,"fastqc")
        args["threads"] = "16"
        _fastqs = fastq_object.get_fastqs()
        not_complete = []
        for fastq in _fastqs:
            sample = os.path.splitext(os.path.splitext(os.path.split(fastq)[1])[0])[0]
            html = os.path.join(args["outdir"],"{}_fastqc.html".format(sample))
            if not os.path.exists(html):
                not_complete.append(fastq)
            else:
                print(html)
        cmds = FastQC.cmd(not_complete,args)
        print(len(not_complete), "FastQC reports not complete...")
        for cmd in cmds:
            print(" ".join(cmd))
            subprocess.call(cmd)
