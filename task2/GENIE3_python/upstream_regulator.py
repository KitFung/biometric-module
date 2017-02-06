import pandas as pd
import scipy.stats as stats

# expression file
fname = '../unknown_expression_data.tsv'


def read_genename(fname):
    with open(fname) as f:
        gene_names = f.readline()
        gene_names = gene_names.strip().split('\t')
    total_gene = len(gene_names)
    return gene_names, total_gene


gene_names, total_gene = read_genename(fname)

# regulatory network

fname = 'VIM.csv'


def read_network(fname, gene_names):
    VIM = pd.read_csv(fname, sep=',', index_col=None, header=None)
    VIM.columns = gene_names
    VIM.index = gene_names
    return VIM


VIM = read_network(fname, gene_names)


# condition(fname):

fname = '../condition.csv'


def read_condition(fname):
    condition = pd.read_csv(fname, sep=',', index_col=0)
    condition.columns = ['cond']
    return condition


condition = read_condition(fname)


# regulate direction

fname = '../regulateDirection.csv'


def read_direction(fname):
    direction = pd.read_csv(fname, sep=',', index_col=0)
    condition.columns = ['cond']
    return condition


direction = read_direction(fname)


# correlation

# predicted state


# overlap p-value.
# p-value < 0.01 = significant

# sub-sample a smaller dataset to simulate
# - left smaller dataset
# - top regulator regulated gene

# use > 0.01 as a threshold to filter out the weak relationship
significant_TF = []

for col in VIM:
    overlap_gene = len(VIM[VIM[col] > 0].index)
    remain_gene = total_gene - overlap_gene

    oddsratio, pvalue = stats.fisher_exact(
        [[1, 1], [remain_gene, overlap_gene]])

    if pvalue < 0.01:
        significant_TF.append((col, pvalue))

print len(gene_names)
print len(significant_TF)

# regulatino direction
# - correlation between 2 gene
# Observed gene regulation
# - expression average different
