import random
import os
import glob
import collections

class SampleSheet(object):

    def __init__(self, filename = None, delim=","):
        if filename:
            attributes = collections.defaultdict(list)
            rows = open(filename,"r").read().splitlines()
            attribute_names = rows.pop(0).split(delim)
            for row in rows:
                for attr, value in zip(attribute_names, row.split(",")):
                    attributes[attr].append(value)
            for attr, value in attributes.items():
                setattr(self, attr.lower(), value)


    def __add__(self, other):
        for attr, value in other.__dict__.items():
            if attr not in self.__dict__:
                setattr(self, attr.lower(), value)
            else:
                self.__dict__[attr] += value
        return self


class FastQDirectory(object):

    def __init__(self, directory):

        self.path = directory
        self.id = "run_{}".format(random.randint(0,100000))
        self.output = os.path.join(self.path, "{}/outs/".format(self.id))
        self.samples = self.get_samples(directory)

    def get_samples(self, directory):
        sheets = glob.glob(os.path.join(directory, "*.csv"))
        sheets += glob.glob(os.path.join(directory, "*.tsv"))
        sample_sheet = SampleSheet()
        for sheet in sheets:
            sample_sheet += SampleSheet(filename=sheet)
        return sample_sheet


    def __eq__(self, other):
        return self.path == self.other

    def out(self):
        return self.output
