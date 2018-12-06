import os
import yaml
yaml_file = os.path.join(os.path.split(os.path.realpath(__file__))[0], "../settings.yaml")

def yaml_configuration():
    with open(yaml_file, "r") as f:
        doc = yaml.load(f)
        return doc

class Configuration(object):
    def __init__(self):
        self.reference = "/igo_large/reference/GRCh38"
        self.rho_matrix = "/home/nceglia/codebase/refdata/rho.json"
        self.copy_number_data = None
        self.scviz_embedding = None
        self.build = self.reference.split("/")[-1]
        self.chemistry = "SC5P-R2"
        overrides = yaml_configuration()
        for attr, value in overrides.items():
            setattr(self, attr, value)
