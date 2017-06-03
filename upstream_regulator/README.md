# Introduction

[Based paper](http://pages.ingenuity.com/rs/ingenuity/images/0812%20upstream_regulator_analysis_whitepaper.pdf)

Try to find the upstream regulator amoung a part of the gene in the dataset using the method provided by `ingenuity.com`.

The process of finding the upstream regulator:

- calculate the overlap p-value of the regulator to the entire dataset
- only use the gene that p-value < 0.01 as regulator
- calculating the gene differentation and regulation
- using those information to fing the activating state of the regulator
- calculate the bias of each regulator
- if the bias is < 0.25
  - calculate the z-score
- else
  - calculate the bias-corrected z-score
- Select the regulator that having absolute value of z-score > 2 as upstream regulator

# Dependency

R library
- DESeq

Python library
- scipy
- pandas
- sklearn
- numpy

# Input:
- count table
  - a matrix
  - have col and row header
  - col: sample, row: gene
- classification of sample
  - a matrix
  - have col and row header
  - clasifed the sample to 2 group ('condA' or 'condB')
  - row: sample, only one row: the condition  ('condA' or 'condB')
- the gene names of entire dataset (optional)
  - using a line break to separate each names

# Output

- the list of upstream regulator (result.csv)
  - list of gene names and its z-score



Example
---

```bash
source pipeline.sh ../example/raw_countTable.csv ../example/condition.csv ../example/result/
```
