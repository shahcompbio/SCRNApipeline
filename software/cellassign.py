import urllib.request
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from urllib.error import HTTPError


class CellAssign(object):

    def __init__(self):
        try:
            bioc_url = urllib.request.urlopen('https://github.com/Irrationone/cellassign/blob/master/R/cellassign_em.R')
            string = ''.join(bioc_url.readlines())
            cellassign_em = SignatureTranslatedAnonymousPackage(string, "cellassign_em")
        except HTTPError:
            pass


    def run_em(self, sce_experiment):
        return None
