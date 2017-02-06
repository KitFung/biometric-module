import csv
from GENIE3 import *

fname = '../unknown_expression_data.tsv'
data = loadtxt(fname, skiprows=1)
f = open(fname)
gene_names = f.readline()
f.close()

gene_names = gene_names.rstrip('\n').split('\t')
(VIM, prediction_score, treeEstimators) = GENIE3(data)

get_link_list(VIM, gene_names=gene_names, file_name='ranking.txt')

with open("VIM.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(VIM)
