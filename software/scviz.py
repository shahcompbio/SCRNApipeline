import subprocess

class SCViz(object):

    @staticmethod
    def cmd(routine, args):
        cmd = ["scvis",routine]
        for flag, value in args.items():
            cmd.append("--{}".format(flag))
            cmd.append(value)
        return cmd

    @staticmethod
    def train(matrix, directory, labels):
        args = dict()
        args["data_matrix_file"] = matrix
        args["out_dir"] = directory
        args["data_label_file"] = labels
        cmd = SCViz.cmd("train", args)
        subprocess.call(cmd)

    @staticmethod
    def map(matrix, directory, embedding):
        args = dict()
        args["data_matrix_file"] = matrix
        args["out_dir"] = directory
        args["pretrained_model_file"] = embedding
        cmd = SCViz.cmd("map",args) 
        subprocess.call(cmd)
