from sklearn.linear_model import LogisticRegression
from scipy.stats.distributions import chi2
from sklearn.metrics import log_loss
from statsmodels.stats import multitest
import operator
import numpy
import itertools
import collections
import tqdm
import os

from software.kallisto import Kallisto
from interface.tenxanalysis import TenxAnalysis
from interface.fastqdirectory import FastQDirectory
from utils.cloud import TenxDataStorage

class DifferentialExpression(object):

    def __init__(self, sampleids, fastqs, chem="v2"):
        assert len(sampleids) == len(fastqs) == 2
        print (sampleids)
        self.model = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
        matrices = dict()
        self.samples = sampleids
        for sampleid, fastq in zip(sampleids, fastqs):
            tenx = TenxDataStorage(sampleid,version="v2")
            tenx.download()
            tenx_analysis = TenxAnalysis(tenx.tenx_path)
            tenx_analysis.load()
            tenx_analysis.extract()
            self.krunner = Kallisto(fastqs, tenx_analysis, chem=chem)
            self.krunner.run_pseudo()
            self.krunner.run_bus()
            matrix = self.krunner.design_matrix()
            matrices[sampleid] = matrix
        self.matrices = matrices
        self.matrix1 = self.matrices[sampleids[0]]
        self.matrix2 = self.matrices[sampleids[1]]
        self.common_genes = set(self.matrix1.keys()).intersection(set(self.matrix2.keys()))


    def logistic_regression(self):
        de_file = "{}_{}_de.tsv".format(self.samples[0],self.samples[1])
        if not os.path.exists(de_file):
            return
            output = open("{}_{}_de.tsv".format(self.samples[0],self.samples[1]),"w")
            output.write("Gene\tPValue\n")
            differential_genes = dict()
            for gene in tqdm.tqdm(self.common_genes):
                tcc_common = set(self.matrix1[gene].keys()).intersection(set(self.matrix2[gene].keys()))
                if len(tcc_common) == 0:
                    continue
                Y = []
                X = []
                cells1 = list(itertools.chain.from_iterable([list(self.matrix1[gene][tcc].keys()) for tcc in tcc_common]))
                cells2 = list(itertools.chain.from_iterable([list(self.matrix2[gene][tcc].keys()) for tcc in tcc_common]))
                if len(cells1) == 0 or len(cells2) == 0:
                    continue
                for cell in cells1:
                    Y.append(self.samples[0])
                    predictors = []
                    for tcc in tcc_common:
                        try:
                            predictors.append(self.matrix1[gene][tcc][cell])
                        except KeyError:
                            predictors.append(0)
                    X.append(predictors)
                for cell in cells2:
                    Y.append(self.samples[1])
                    predictors = []
                    for tcc in tcc_common:
                        try:
                            predictors.append(self.matrix2[gene][tcc][cell])
                        except KeyError:
                            predictors.append(0)
                    X.append(predictors)
                classes = set(Y)
                Y = numpy.array(Y)
                X = numpy.array(X)
                if Y.shape[0] < 2 or len(classes) == 1:
                    continue
                self.model.fit(X,Y)
                null_prob = 2.0 / float(Y.shape[0]) * numpy.ones(Y.shape)
                df = X.shape[1]
                alt_prob = self.model.predict_proba(X)
                alt_log_likelihood = -log_loss(Y, alt_prob, normalize=False)
                null_log_likelihood = -log_loss(Y, null_prob, normalize=False)
                G = 2 * (alt_log_likelihood - null_log_likelihood)
                p_value = chi2.sf(G, df)
                differential_genes[gene] = p_value
                output.write("{}\t{}\n".format(gene,p_value))
            sorted_genes = sorted(differential_genes.items(), key=operator.itemgetter(1))
            print ("**************** Differential Genes ********************")
            for gene, pvalue in sorted_genes[:100]:
                print (gene, pvalue)
            output.close()
        else:
            differential_genes = dict()
            differential_genes_adj = dict()
            genes = open(de_file,"r").read().splitlines()
            genes.pop(0)
            _genes = []
            pvalues = []
            adjpvalues = []
            for gene in genes:
                gene, pvalue = gene.split()
                differential_genes[gene] = float(pvalue)
                pvalues.append(float(pvalue))
                _genes.append(gene)
            adj_pvalues = list(multitest.multipletests(pvalues)[1])
            print(adj_pvalues)
            for gene, pvalue, adjp in zip(_genes,pvalues, adj_pvalues):
                differential_genes_adj[gene] = adjp
            sorted_genes = sorted(differential_genes_adj.items(), key=operator.itemgetter(1))
            thresholds = (0.05,0.01,0.001)
            import collections
            sig_genes = collections.defaultdict(list)
            for gene, pvalue in sorted_genes:
                for threshold in thresholds:
                    if pvalue < threshold:
                        sig_genes[str(threshold)].append(gene)
            print ("**************** Differential Genes ********************")
            for thresh, sig_genes in sig_genes.items():
                print (thresh, len(sig_genes))
            for gene, pvalue in sorted_genes[:100]:
                print (gene, pvalue)
        return sorted_genes


if __name__ == '__main__':
    from itertools import combinations
    samples = ["Y7640", "Y7652", "Y7668", "Y8841"]
    fastqs = {}
    for sample in samples:
        output = "/igo_large/scratch/de_ciara/"
        fastq_directory = FastQDirectory("/igo_large/data/{}".format(sample), sample, output)
        fastqs[sample] = fastq_directory

    pairwise = combinations(samples,2)
    for pair in pairwise:
        fastq_set = [fastqs[pair[0]],fastqs[pair[1]]]
        de = DifferentialExpression(pair, fastq_set)
        res = de.logistic_regression()
