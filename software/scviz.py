from interface.singlecellexperiment import SingleCellExperiment
import pickle
import subprocess

class SCViz(object):

    @staticmethod
    def cmd(routine, args):
        cmd = ["scvis",routine]
        for flag, value in args.items():
            cmd.append("--{}".format(flag))
            cmd.append(value)
        cmd.append("--verbose")
        return cmd

    @staticmethod #TODO move to TenxAnalysis
    def create_input_files(sce_experiment, tenxanalysis, cellassignments):
        labels = tenxanalysis.clustering_labels()
        barcodes = sce_experiment.colData["Barcode"]
        # assert len(barcodes) == len(cellassignments), "Provided cell assignments do not cover all barcodes."
        # assert len(set(barcodes).difference(set(labels.keys()))) == 0, "Not all barcodes have cluster labels."
        assignments = dict(zip(barcodes,cellassignments))
        output = open("labels.tsv","w")
        output.write("Cluster\n")
        for barcode, cluster in labels.items():
            try:
                output.write(assignments[barcode] + "\n")
            except KeyError as e:
                continue
        output.close()
        tenxanalysis.map_projection(assignments, "matrix.tsv")

    @staticmethod
    def train(rdata, tenxanalysis, cellassignments):
        print("Loading SCVIS")
        fit = pickle.load(open(cellassignments,"rb"))
        pyfit = dict(zip(fit.names, list(fit)))
        cellassignments = pyfit["cell_type"]
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        SCViz.create_input_files(sce_experiment, tenxanalysis, cellassignments)
        args = dict()
        args["data_matrix_file"] = "matrix.tsv"
        args["out_dir"] = "./"
        args["data_label_file"] = "labels.tsv"
        #args["normalize"] = "1.0"
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
