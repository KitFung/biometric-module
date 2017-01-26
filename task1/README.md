# Introduction

Doing consensus cluster and classification to a gene expression matrix file - each row represent a gene, each column represent a sample, each field represent the expression value of that gene in that sample.

*The processing step of clustering*

1. Filtering the gene (only use the gene that have no zero in all sample)
2. Using `YuGene` to remove batch effect and normalization
3. Compute the gap statistic to find the best K (try 1 to 6 K)
4. Using median absolute deviation(>0.01) to find the most valuable genes
5. Do consensus clustering base on it and calculate silhouette

*The processing step of classification*

1. Filtering the gene (only use the gene that have no zero in all sample)
2. Using `YuGene` to remove batch effect and normalization both train and test dataset
3. Finding the common gene between training and testing dataset. (The training and testing dataset must have at least 2 common gene)
4. Search for top differential genes between each subtypes and use it in later classify (If the number of gene after finding is least than 5, just skip this)
5. Do the filtering to the gene with condtion (AUC > 0.9) and use it in later classify (If the number of gene after filtering is least than 5, just skip this)
6.  To select the optimal threshold for centroid shrinkage, we performed cross-validation to select the one yielding a good performance with the least number of genes
7. Create a classifier([nearest shrunken centroid classification](http://statweb.stanford.edu/~tibs/PAM/Rdist/howwork.html)) using the parameter found in previous step
8. Classify the test data using the classifier

[The related experiment](https://bioconductor.org/packages/release/data/experiment/vignettes/DeSousa2013/inst/doc/DeSousa2013-Vignette.pdf)
[The related paper](http://dare.uva.nl/document/2/127041)

# Dependency

### R Library

- [YuGene](https://cran.r-project.org/web/packages/YuGene/index.html)

- [BioBase](https://bioconductor.org/packages/release/bioc/html/Biobase.html)

- [cluster](https://cran.r-project.org/web/packages/cluster/cluster.pdf)

- [siggenes](http://bioconductor.org/packages/release/bioc/html/siggenes.html)

- [ConsensusClusterPlus](https://www.bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html)

- [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)

- [pamr](https://cran.r-project.org/web/packages/pamr/index.html)

- [gplots](https://cran.r-project.org/web/packages/gplots/index.html)

- [ROCR](https://cran.r-project.org/web/packages/ROCR/ROCR.pdf)

# Usage

```bash
mkdir output_folder
# clustering
Rscript clustering.R gene_data_matrix.txt output_folder
# classification
Rscript classification.R training_dataset.csv training_label.csv test_dataset.csv output_folder
```

# Input Data Format

The columns is sample, the row is gene.
The expression data is separated by tab
The label file is in format `"gene",class_number` with header

# Output File

### In clustering

- **`consesClusterResult.csv` - the clustering result**
- `gap.pdf` - the "gap" value while using different K cluster, it measure goodness of clustering
- `silh.pdf` - the silhouette value while using different K cluster. `Silhouette refers to a method of interpretation and validation of consistency within clusters of data`
- `icl.pdf` - the result after calculating cluster-consensus and item-consensus.

### In classification

- **`classifed_result.csv`** - the classification result
- `diffGenes.f.csv` - the gene used in the classifier
- `thresh.csv` - the threshold value obtained from step 6
- `centroids.csv` - A matrix of (unshrunken) class centroids, n by nclass
- ``centroid.overall.csv` - A vector containing the (unshrunken) overall centroid (all classes together)
