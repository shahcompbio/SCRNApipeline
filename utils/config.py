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
jobpath: "{2}"
datapath: "{3}"
referencepath: "{4}"
rho_matrix: /work/shah/reference/transcriptomes/markers/hgsc_v1.yaml
cellranger: "{5}"
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
qc_type: "standard"
scviz_embedding: None
copy_number_data: None
mito: 10
lsf: {6}
"""


def write_config(prefix, build, jobpath, datapath, referencepath, cellranger, lsf):
    output = open("settings.yaml","w")
    output.write(basic_yaml.format(prefix, build, jobpath, datapath, referencepath, cellranger, lsf))

class Configuration(object):
    def __init__(self):
        self.build = "GRCh38"
        overrides = yaml_configuration()
        if overrides != None:
            for attr, value in overrides.items():
                setattr(self, attr, value)
        refobj = ReferenceDataStorage(self.build, self.referencepath)
        if not hasattr(self, "reference"):
            self.reference = refobj.download()
        self.genes_gtf = os.path.join(self.reference, "genes/genes.gtf")
