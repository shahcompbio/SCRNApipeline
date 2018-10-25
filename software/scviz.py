from interface.singlecellexperiment import SingleCellExperiment
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
    def create_input_files(sce_experiment, tenxanalysis, cellassignments):
        labels = tenxanalysis.clustering_labels(graphclust=True)
        barcodes = sce_experiment.colData["Barcode"]
        assert len(barcodes) == len(cellassignments), "Provided cell assignments do not cover all barcodes."
        assert set(labels.keys()) == set(barcodes), "Not all barcodes have cluster labels."
        assignments = dict(zip(barcodes,cellassignments))
        output = open("datamatrix.tsv","w")
        output.write("Cluster\n")
        for barcode, cluster in labels.items():
            output.write(cluster + "\n")
        output.close()
        tenxanalysis.map_projection(assignments, "labels.tsv")

    @staticmethod
    def train(rdata, tenxanalysis, cellassignments):
        fit = pickle.load(open("cell_assign_fit.pkl","rb"))
        pyfit = dict(zip(fit.names, list(fit)))
        cellassignments = pyfit["cell_type"]
        sce_experiment = SingleCellExperiment.fromRData(rdata)
        SCViz.create_input_files(sce_experiment, tenxanalysis, cellassignments)
        args = dict()
        args["data_matrix_file"] = "datamatrix.tsv"
        args["out_dir"] = "./"
        args["data_label_file"] = "labels.tsv"
        cmd = SCViz.cmd("train", args)
        print("Running...")
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
