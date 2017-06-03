import sys
import os
import math
import csv
import pandas as pd
import scipy.stats as stats


def read_exprs(fname):
    exprs = pd.read_csv(fname, sep='\t', index_col=None)
    gene_names = list(exprs.columns)
    total_gene = len(gene_names)
    return exprs, gene_names, total_gene


def read_network(fname, gene_names):
    VIM = pd.read_csv(fname, sep=',', index_col=None, header=None)
    VIM.columns = gene_names
    VIM.index = gene_names
    return VIM


def read_observed_regulation(fname):
    observed_regulation = pd.read_csv(fname, sep=',', index_col=0)
    observed_regulation.columns = ['cond']
    return observed_regulation


def read_names(fname):
    names = []
    with open(fname) as f:
        names = [line.strip() for line in f.readlines()]
    return names


# correlation -> regulate direction
# > 0 = activating
# <= 0 = inhibiting
def cal_direction(exprs):
    correlation = exprs.corr()
    direction = correlation
    direction[direction > 0] = 1
    direction[direction <= 0] = -1
    return direction


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


def divide_TR_by_bias(TRs, bias):
    unbiased_TRs = []
    biased_TRs = []
    for TR in TRs:
        if bias[TR] < 0.25:
            unbiased_TRs.append(TR)
        else:
            biased_TRs.append(TR)
    return unbiased_TRs, biased_TRs


# z-score
def cal_zscore(TRs, TR2TGs, weight, predicted_state):
    z_scores = []
    for TR in TRs:
        z_score = 0
        TGs = TR2TGs[TR]
        square_sum = 0
        for TG in TGs:
            z_score += weight.loc[TR, TG] * predicted_state.loc[TR, TG]
            square_sum += weight.loc[TR, TG] ** 2
        z_score /= math.sqrt(square_sum)
        z_scores.append((TR, z_score))
    z_scores.sort(key=lambda t: t[1], reverse=True)
    return z_scores


# bias-corrected z-score
def cal_bias_corrected_zscore(TRs, TR2TGs, weight, predicted_state, bias):
    z_scores = []
    for TR in TRs:
        z_score = 0
        TGs = TR2TGs[TR]
        square_sum = 0
        for TG in TGs:
            tmp = predicted_state.loc[TR, TG] - bias[TR]
            z_score += weight.loc[TR, TG] * tmp
            square_sum += weight.loc[TR, TG] ** 2
        z_score /= math.sqrt(square_sum)
        z_scores.append((TR, z_score))
    z_scores.sort(key=lambda t: t[1], reverse=True)
    return z_scores


if __name__ == '__main__':
    exprs_fname = None  # '../unknown_expression_data.tsv'
    VIM_fname = None  # 'VIM.csv'
    observed_regulation_fname = None  # '../observedRegulation.csv'
    save_folder = None
    allgene_fname = None  #
    sub_sampling_genes = None
    if len(sys.argv) == 5:
        exprs_fname = sys.argv[1]
        VIM_fname = sys.argv[2]
        observed_regulation_fname = sys.argv[3]
        save_folder = sys.argv[4]
    elif len(sys.argv) == 6:
        exprs_fname = sys.argv[1]
        VIM_fname = sys.argv[2]
        observed_regulation_fname = sys.argv[3]
        save_folder = sys.argv[4]
        allgene_fname = sys.argv[5]
    else:
        print("""
              Usage: python %s
              exprs_fname
              VIM_fname
              observed_regulation_fname
              allgene_fname
              """ % sys.argv[0])
        exit()

    # expression
    exprs, gene_names, total_gene = read_exprs(exprs_fname)
    # regulatory network
    VIM = read_network(VIM_fname, gene_names)
    # Observed gene regulation
    observed_regulation = read_observed_regulation(observed_regulation_fname)

    direction = cal_direction(exprs)
    predicted_state = cal_predicted_state(VIM, direction, observed_regulation)

    TRs = VIM.index

    if allgene_fname is None:
        dataset_genes = gene_names[::2]
    else:
        dataset_genes = read_names(allgene_fname)

    TR2TGs = form_TR2TGs(TRs, VIM)
    # filtering using p-value < 0.01
    TRs_after_pvalue = filter_TR_by_pvalue(
        TRs, TR2TGs, VIM, dataset_genes)

    bias_data = cal_bias_data(observed_regulation)
    bias = cal_bias(TRs_after_pvalue, TR2TGs,  direction, bias_data)

    unbiased_TRs, biased_TRs = divide_TR_by_bias(TRs_after_pvalue, bias)

    # if bias < 0.25, use cal_zscore else use bias_corrected_zscore
    z_scores = cal_zscore(unbiased_TRs, TR2TGs, VIM, predicted_state)
    bias_corrected_zscore = cal_bias_corrected_zscore(
        biased_TRs, TR2TGs, VIM, predicted_state, bias)

    # combine result
    z_scores.extend(bias_corrected_zscore)
    z_scores.sort(key=lambda t: t[1], reverse=True)

    save_fname = os.path.join(save_folder, 'result.csv')
    significant_result = [(gene, val)
                          for gene, val in z_scores if abs(val) > 2]
    with open(save_fname, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(significant_result)
