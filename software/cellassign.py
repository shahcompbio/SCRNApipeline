import urllib.request
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage


bioc_url = urllib.request.urlopen('https://github.com/Irrationone/cellassign/blob/master/R/cellassign_em.R')
string = ''.join(bioc_url.readlines())

cellassign_em = SignatureTranslatedAnonymousPackage(string, "cellassign_em")
