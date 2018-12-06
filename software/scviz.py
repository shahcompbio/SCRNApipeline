from interface.singlecellexperiment import SingleCellExperiment
from interface.tenxanalysis import TenxAnalysis
import scanpy.api as pl
import pickle
import subprocess
import tqdm
import numpy
import os
import scipy

class SCViz(object):

    @staticmethod
    def cmd(routine, args):
        cmd = ["scvis",routine]
        for flag, value in args.items():
            cmd.append("--{}".format(flag))
            cmd.append(value)
        cmd.append("--verbose")
        print(" ".join(cmd))
        return cmd

    @staticmethod
    def generate_config(perplexity, components, output):
        filename = os.path.join(output,"model.yaml")
        yaml = SCViz.config(perplexity, components)
        output = open(filename,"w")
        output.write(yaml)
        output.close()
        return filename

    @staticmethod
    def create_input_files(rdata, components, output):
        sce = SingleCellExperiment.fromRData(rdata)
        embedding = sce.reducedDims["PCA"]
        counts = []
        for i in range(0,len(embedding),components):
            counts.append(embedding[i:i+(components)])
        print(counts)
        for row in counts:
            print(len(row))
        print(len(counts))
        counts = numpy.array(counts)
        print(counts.shape)
        header = []
        for c in range(components):
            header.append("PC_{}".format(c))
        header = "\t".join(header)
        filename = os.path.join(output,"matrix.tsv")
        numpy.savetxt(filename, counts, delimiter="\t", header=header)
        return filename

    @staticmethod
    def train(output, perplexity, components, rdata):
        output = os.path.join(output, "{}_{}_raw".format(perplexity, components))
        try:
            os.makedirs(output)
        except Exception as e:
            pass
        yaml = SCViz.generate_config(perplexity, components, output)
        matrix_file = SCViz.create_input_files(rdata, components, output)
        args = dict()
        args["data_matrix_file"] = matrix_file
        args["out_dir"] = output
        args["config_file"] = yaml
        cmd = SCViz.cmd("train", args)
        print("Running...\n")
        print(" ".join(cmd)+"\n")
        subprocess.call(cmd)
        print("Finished!")

    @staticmethod
    def map(matrix, embedding, output):
        args = dict()
        args["data_matrix_file"] = matrix
        args["out_dir"] = directory
        args["pretrained_model_file"] = embedding
        cmd = SCViz.cmd("map",args)
        subprocess.call(cmd)

    @staticmethod
    def config(perplexity, components):
        configuration = """
        hyperparameter: {
          optimization: {
            method: Adam,
            learning_rate: 0.01
          },

          batch_size: 512,
          max_epoch: 50,
          regularizer_l2: 0.001,

          perplexity: """
        configuration += str(perplexity) + ", "
        configuration += """
          seed: 1
        }

        architecture: {
          latent_dimension: """
        configuration +=  "2, "
        configuration += """
          inference: {
            layer_size: [128, 64, 32],
          },

          model: {
            layer_size: [32, 32, 32, 64, 128],
          },

          activation: "ELU"
        }
        """
        return configuration


if __name__ == '__main__':
    valid_barcodes = list(SingleCellExperiment.fromRData("/home/nceglia/jobs/post_treatment_full/sce_final.rdata").colData["Barcode"])
    valid_symbols = list(SingleCellExperiment.fromRData("/home/nceglia/jobs/post_treatment_full/sce_final.rdata").rowData["Symbol"])
    analysis = TenxAnalysis("/home/nceglia/jobs/post_treatment_full/run_post_treatment_full/outs")
    pkl_fit = "/home/nceglia/jobs/post_treatment_full/cell_assign_fit.pkl"
    fit = pickle.load(open(pkl_fit,"rb"))
    SCViz.create_input_files(analysis, fit, None, valid_barcodes, valid_symbols)
