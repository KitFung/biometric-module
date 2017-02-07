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
def form_TR2TGs(TRs, VIM):
    TR2TGs = {}
    for TR in TRs:
        TGs_case = VIM.loc[TR] > 0
        TGs = TGs_case[TGs_case].index
        TR2TGs[TR] = TGs
    return TR2TGs


def cal_overlap_pvalue(all_genes, regulated_gene):
    all_genes = set(all_genes)
    regulated_gene = set(regulated_gene)
    oddsratio, pvalue = stats.fisher_exact(
        [[0,                                len(regulated_gene - all_genes)],
         [len(all_genes - regulated_gene),  len(all_genes & regulated_gene)]])
    return pvalue


def filter_TR_by_pvalue(TRs, TR2TGs, VIM, dataset_genes):
    significant_TRs = []
    for TR in TRs:
        TGs = TR2TGs[TR]
        pvalue = cal_overlap_pvalue(dataset_genes, TGs)
        if pvalue < 0.01:
            significant_TRs.append(TR)
    return significant_TRs


TRs = VIM.index
sub_sampling_genes = gene_names[::2]
TR2TGs = form_TR2TGs(TRs, VIM)
TRs_after_pvalue = filter_TR_by_pvalue(TRs, TR2TGs, VIM, sub_sampling_genes)


# z-score
def cal_zscore(TRs, TR2TGs, weight, predicted_state):
    z_scores = {}
    for TR in TRs:
        z_scores[TR] = 0
        TGs = TR2TGs[TR]
        for TG in TGs:
            z_scores[TR] += weight.loc[TR, TG] * predicted_state.loc[TR, TG]
    return z_scores


z_scores = cal_zscore(TRs, TR2TGs, VIM, predicted_state)


# bias = bias_data * bias_TR
def cal_bias_data(observed_regulation):
    N_up = observed_regulation[observed_regulation > 0].count()['cond']
    N_down = observed_regulation[observed_regulation < 0].count()['cond']
    return float(N_up - N_down)/(N_up + N_down)


def cal_bias_TR(TR, TGs, direction):
    tTR = direction.loc[TR, TGs]
    activating = tTR[tTR > 0].count()
    inhibiting = tTR[tTR < 0].count()
    return float(activating - inhibiting)/(activating + inhibiting)


def cal_bias(TRs, TR2TGs, direction, bias_data):
    bias = {}
    for TR in TRs:
        bias[TR] = bias_data * cal_bias_TR(TR, TR2TGs[TR], direction)
    return bias


def filter_TR_by_bias(TRs, bias):
    significant_TRs = []
    for TR in TRs:
        if bias[TR] < 0.25:
            significant_TRs.append(TR)
    return significant_TRs


bias_data = cal_bias_data(observed_regulation)
bias = cal_bias(TRs_after_pvalue, TR2TGs,  direction, bias_data)
TRs_after_bias = filter_TR_by_bias(TRs_after_pvalue, bias)
