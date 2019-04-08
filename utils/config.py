import os
import yaml

from utils.cloud import ReferenceDataStorage

yaml_file = os.path.join(os.getcwd(), "settings.yaml")
def yaml_configuration():
    if os.path.exists(yaml_file):
        with open(yaml_file, "r") as f:
            doc = yaml.load(f)
            return doc


basic_yaml = """
prefix: "{0}"
build: "{1}"
jobpath: "/jobs/test"
datapath: "/data/test"
reference: "{2}"
rho_matrix: None
cellranger: "{3}"
copy_number_data: None
scviz_embedding: None
run_scvis: False
run_cellassign: False
run_clonealign: False
plot_scvis: False
clustering: False
report: True
perplexity: 5
resolution: 0.2
stds: 6
components: 50
chemistry: auto
low_counts_genes_threshold: 4
"""

cellranger = "/home/ceglian/codebase/cellranger-2.2.0/cellranger-cs/2.2.0/bin"

def write_config(prefix, output, reference):
    output = open("settings.yaml","w")
    if "/" not in reference:
        build = reference
        reference = None
    else:
        build = reference.split("/")[-2]
        reference = reference
    output.write(basic_yaml.format(prefix, build, reference, cellranger))

class Configuration(object):
    def __init__(self):
        self.reference = None
        self.copy_number_data = None
        self.scviz_embedding = None
        self.datapath = None
        self.jobpath = None
        self.referencepath = None
        self.build = None
        self.chemistry = "auto"
        self.qc_type = "standard"
        self.low_counts_genes_threshold = 4
        self.stds = 6
        overrides = yaml_configuration()
        if overrides != None:
            for attr, value in overrides.items():
                setattr(self, attr, value)
        if self.reference == "None" or self.reference == None:
            if self.referencepath == None:
                self.referencepath = "/reference"
            refobj = ReferenceDataStorage(self.build, self.referencepath)
            self.reference = refobj.download()
        else:
            self.genes_gtf = os.path.join(self.reference, "genes/genes.gtf")
        if self.build == None:
            self.build = "GRCh38"
