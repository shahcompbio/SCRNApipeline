from interface.singlecellexperiment import SingleCellExperiment
import pickle
import subprocess
import numpy

class SCViz(object):

    @staticmethod
    def cmd(routine, args):
        cmd = ["scvis",routine]
        for flag, value in args.items():
            cmd.append("--{}".format(flag))
            cmd.append(value)
        cmd.append("--verbose")
        return cmd

    @staticmethod
    def create_input_files(sce_experiment, cell_assign_fit):
        cell_types = dict(zip(cell_assign_fit["Barcode"], cell_assign_fit["cell_type"]))
        barcodes = sce_experiment.colData["Barcode"]
        counts = sce_experiment.assays["counts"].transpose()
        data_matrix = []
        output = open("labels.tsv","w")
        output.write("CellType\n")
        for barcode, features in zip(barcodes, counts.toarray()):
            try:
                cell_type = cell_types[barcode]
                output.write(cell_type + "\n")
                data_matrix.append(list(features))
            except KeyError as e:
                continue
        output.close()
        header = "\t".join(sce_experiment.rowData["Symbol"])
        matrix = numpy.array(data_matrix)
        print(matrix.shape)
        assert len(sce_experiment.rowData["Symbol"]) == matrix.shape[1]
        numpy.savetxt("matrix.tsv",matrix,delimiter="\t",header=header)

    @staticmethod
    def train(rdata, output, cellassignments):
        print("Loading SCVIS")
        fit = pickle.load(open(cellassignments,"rb"))
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        SCViz.create_input_files(sce_experiment, fit)
        args = dict()
        args["data_matrix_file"] = "matrix.tsv"
        args["out_dir"] = output
        args["data_label_file"] = "labels.tsv"
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
