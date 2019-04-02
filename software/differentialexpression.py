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

import scanpy.api as sc
from software.kallisto import Kallisto
from interface.tenxanalysis import TenxAnalysis
from interface.fastqdirectory import FastQDirectory
from software.batchcorrection import Scanorama
from utils.cloud import TenxDataStorage
from interface.qualitycontrol import QualityControl

class DifferentialExpression(object):

    def __init__(self, sampleids, chem="v2", output="./"):
        self.output = output
        self.samples = sampleids
        self.tenxs = []
        for sampleid in self.samples:
            tenx = TenxDataStorage(sampleid,version=chem)
            tenx.download()
            tenx_analysis = TenxAnalysis(tenx.tenx_path)
            tenx_analysis.load()
            tenx_analysis.extract()
            self.tenxs.append(tenx_analysis)

    def run(self, method="logreg", correction="benjamini-hochberg", n_genes=500, mito=0.1, min_genes=200, min_cells=3):
        genesets = []
        adatas = []
        for sampleid, tenx in zip(self.samples, self.tenxs):
            adata = tenx.create_scanpy_adata_basic(sample_key=self.samples[0])
            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            sc.pp.filter_cells(adata, min_genes=200, inplace=True)
            sc.pp.filter_genes(adata, min_cells=3,inplace=True)
            mito_genes = adata.var_names.str.startswith('MT-')
            adata.obs['percent_mito'] = numpy.sum(adata[:, mito_genes].X, axis=1).A1 / numpy.sum(adata.X, axis=1).A1
            adata.obs['n_counts'] = adata.X.sum(axis=1).A1
            adata = adata[adata.obs['percent_mito'] < mito, :]
            sc.pp.normalize_per_cell(adata)
            sc.pp.log1p(adata)
            adata.raw = adata
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata = adata[:, adata.var['highly_variable']]
            sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
            sc.pp.scale(adata, max_value=10)
            common_genes = list(adata.var.index)
            genesets.append(set(common_genes))
            adatas.append(adata)
        _adatas = []
        common_genes = set.intersection(*genesets)
        common_genes = list(common_genes)
        for adata in adatas:
            _adatas.append(adata[:,common_genes])
        corrected = Scanorama.integrate_and_correct(_adatas)
        batch_to_sample  = dict(zip(range(len(self.samples)),self.samples))
        corrected.obs['sample'] = [batch_to_sample[int(i)] for i in corrected.obs['batch']]
        sc.pp.scale(adata, max_value=10)
        sc.tl.rank_genes_groups(corrected, groupby="sample", method='t-test', n_genes=n_genes)
        foldchange = dict()
        pvalues = dict()
        adjpval = dict()
        symbols = dict()
        for i, sample in enumerate(self.samples):
            print ("Writing {}".format(sample))
            matrix = []
            columns = []
            for key, values in corrected.uns['rank_genes_groups'].items():
                if key == "params": continue
                columns.append(key)
                scores = [x[i] for x in list(values)]
                if key == "pvals":
                    pvalues[sample] = list(map(float, scores))
                if key == "pvals_adj":
                    adjpval[sample] = list(map(float, scores))
                if key == "logfoldchanges":
                    foldchange[sample] = list(map(float, scores))
                if key == "names":
                    symbols[sample] = scores
                matrix.append(scores)
            output = open(os.path.join(self.output, "{}_{}_de.tsv".format(method, sample)),"w")
            output.write("\t".join(columns)+"\n")
            for row in zip(*matrix):
                row = map(str,row)
                output.write("\t".join(row)+"\n")
            output.close()
        top_genes = []
        brackets = []
        start = 0
        stop = 10
        stride = 10

        for i, sample in enumerate(self.samples):
            print ("Sorting {}".format(sample))

            pvalue_dict = dict(zip(symbols[sample], pvalues[sample]))
            adj_pvalue_dict = dict(zip(symbols[sample], adjpval[sample]))
            fold_change_dict = dict(zip(symbols[sample], foldchange[sample]))

            sort_dict = fold_change_dict

            filtered_stats = dict()
            for symbol in symbols[sample]:
                if pvalue_dict[symbol] > 0.05: continue
                if adj_pvalue_dict[symbol] > 0.2: continue
                if abs(fold_change_dict[symbol]) < 1.0: continue
                filtered_stats[symbol] = sort_dict[symbol]

            print ("Filtered Symbols",len(filtered_stats.keys()))
            sorted_symbols = list(reversed(sorted(filtered_stats.items(), key=operator.itemgetter(1))))

            top_genes += [x[0] for x in sorted_symbols[start:stop]]

            print("Symbols", sample, [x[0] for x in sorted_symbols[start:stop]])
            print("Fold Changes",sample, [fold_change_dict[x[0]] for x in sorted_symbols[start:stop]])
            print("P-Values",sample, [pvalue_dict[x[0]] for x in sorted_symbols[start:stop]])
            print("Adj P-Values",sample, [adj_pval_dict[sample][x[0]] for x in sorted_symbols[start:stop]])

            brackets.append((start,stop))
            start += stride
            stop += stride

        sc.pl.stacked_violin(corrected, scale='width', var_names=top_genes, groupby="sample", var_group_positions=brackets, var_group_labels=self.samples, save="de.png")

    def run_transcript(self, fastqs=[]):
        matrices = dict()
        assert len(fastqs) == len(self.samples), "Provide fastq object for each sample."
        for sampleid, fastq in zip(self.samples, self.fastqs):
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
        self.model = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
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
