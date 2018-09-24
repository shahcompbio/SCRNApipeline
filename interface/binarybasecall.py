import random

class BinaryBaseCall(object):

    def __init__(self, directory):
        self.path = directory
        self.id = "run_{}".format(random.randint(0,100000))
        self.csv = "/home/ceglian/data/cellranger-tiny-bcl-simple-1.2.0.csv"

    def __eq__(self, other):
        return self.path == other.path
