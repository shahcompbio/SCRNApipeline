import random
import os

class FastQDirectory(object):

    def __init__(self, directory):
        self.path = directory
        self.id = "run_{}".format(random.randint(0,100000))
        self.output = os.path.join(self.path, "{}/outs/".format(self.id))

    def __eq__(self, other):
        return self.path == self.other

    def out(self):
        return self.output
