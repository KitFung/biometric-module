import pandas as pd
import scipy.stats as stats

# expression file
fname = '../unknown_expression_data.tsv'


def read_exprs(fname):
    exprs = pd.read_csv(fname, sep='\t', index_col=None)
    gene_names = list(exprs.columns)
    total_gene = len(gene_names)
    return exprs, gene_names, total_gene


exprs, gene_names, total_gene = read_exprs(fname)


# regulatory network
fname = 'VIM.csv'


def read_network(fname, gene_names):
    VIM = pd.read_csv(fname, sep=',', index_col=None, header=None)
    VIM.columns = gene_names
    VIM.index = gene_names
    return VIM


VIM = read_network(fname, gene_names)


# # condition(fname):
# fname = '../condition.csv'
#
#
# def read_condition(fname):
#     condition = pd.read_csv(fname, sep=',', index_col=0)
#     condition.columns = ['cond']
#     return condition
#
#
# condition = read_condition(fname)


# Observed gene regulation
# 1 = up
# -1 = down
# NA = skip
fname = '../observedRegulation.csv'


def read_observed_regulation(fname):
    observed_regulation = pd.read_csv(fname, sep=',', index_col=0)
    observed_regulation.columns = ['cond']
    return observed_regulation


observed_regulation = read_observed_regulation(fname)


# correlation -> regulate direction
# > 0 = activating
# <= 0 = inhibiting
def cal_direction(exprs):
    correlation = exprs.corr()
    direction = correlation
    direction[direction > 0] = 1
    direction[direction <= 0] = -1
    return direction


direction = cal_direction(exprs)


# predicted state

def cal_predicted_state(VIM, direction, observed_regulation):
    predicted_state = VIM.copy()
    for TR in gene_names:
        for TG in gene_names:
            if VIM.loc[TR, TG] > 0:
                if direction.loc[TR, TG] == \
                   observed_regulation.loc[TG, 'cond']:
                    predicted_state.loc[TR, TG] = 1
                else:
                    predicted_state.loc[TR, TG] = -1
    return predicted_state


predicted_state = cal_predicted_state(VIM, direction, observed_regulation)


# overlap p-value.
# p-value < 0.01 = significant

# sub-sample a smaller dataset to simulate
# - left smaller dataset
# - top regulator regulated gene

# use > 0.01 as a threshold to filter out the weak relationship
def cal_overlap_pvalue(all_genes, regulated_gene):
    all_genes = set(all_genes)
    regulated_gene = set(regulated_gene)
    oddsratio, pvalue = stats.fisher_exact(
        [[0,                                len(regulated_gene - all_genes)],
         [len(all_genes - regulated_gene),  len(all_genes & regulated_gene)]])
    return pvalue


def filter_TR_by_pvalue(VIM, dataset_genes):
    significant_TRs = []
    for TR in VIM.index:
        TG_case = VIM.loc[TR] > 0
        TG = TG_case[TG_case].index
        pvalue = cal_overlap_pvalue(dataset_genes, TG)
        if pvalue < 0.01:
            significant_TRs.append(TR)
    return significant_TRs


sub_sampling_genes = gene_names[::2]
TRs_after_pvalue = filter_TR_by_pvalue(VIM, sub_sampling_genes)


# z-score
def cal_zscore(VIM, TRs, weight, predicted_state):
    for TR in TRs:
        pass
    pass
