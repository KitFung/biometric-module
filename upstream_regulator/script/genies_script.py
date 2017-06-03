import csv
import sys
import os
from GENIE3 import *

fname = sys.argv[1]  # '../unknown_expression_data.tsv'
save_folder = sys.argv[2]

data = loadtxt(fname, skiprows=1)
f = open(fname)
gene_names = f.readline()
f.close()

gene_names = gene_names.rstrip('\n').split('\t')
(VIM, prediction_score, treeEstimators) = GENIE3(data)

with open(os.path.join(save_folder, "VIM.csv"), "wb") as f:
    writer = csv.writer(f)
    writer.writerows(VIM)
