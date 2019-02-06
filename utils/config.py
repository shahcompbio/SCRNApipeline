import os
import yaml
# if os.path.exists("")
# yaml_file = os.path.join(os.path.split(os.path.realpath(__file__))[0], "../settings.yaml")
yaml_file = os.path.join(os.getcwd(), "settings.yaml")
def yaml_configuration():
    with open(yaml_file, "r") as f:
        doc = yaml.load(f)
        return doc

class Configuration(object):
    def __init__(self):
        self.reference = "/igo_large/reference/GRCh38"
        self.copy_number_data = None
        self.scviz_embedding = None
        self.build = self.reference.split("/")[-1]
        self.chemistry = "SC5P-R2"
        self.qc_type = "tenx_filtered"
        overrides = yaml_configuration()
        for attr, value in overrides.items():
            setattr(self, attr, value)
        self.genes_gtf = os.path.join(self.reference, "genes/genes.gtf")
