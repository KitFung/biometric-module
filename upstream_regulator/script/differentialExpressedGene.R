args = commandArgs(TRUE)
COUNT_PATH = args[1]  # 'raw_countTable.csv'
COND_PATH = args[2]  # 'raw_meta.csv'
SAVE_DIR = args[3]  # '../example


# processing raw data
countTable = read.csv(COUNT_PATH, header=TRUE, row.names = 1, check.names=FALSE)
countTable = countTable[rowSums(countTable) != 0, ]
condition = read.csv(COND_PATH, header=TRUE, row.names = 1, check.names=FALSE)
condition = as.vector(condition[, colnames(condition)[1]])

library(DESeq)

# Use DESeq to discover differential expressed geene
# pval : p value for the statistical significance of this change
# padj : p value adjusted for multiple testing with the Benjamini-Hochberg procedure (see the R function
#        p.adjust), which controls false discovery rate (FDR)
# use the padj < 0.1 as condition, it it filtered all gene, use pval < 0.1
differentialExpressedGene <- function(countTable, condition) {

  cds <- newCountDataSet( countTable, condition )
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds )
  res <- nbinomTest( cds, "condA", "condB" )
  rownames(res) <- rownames(countTable)
  resSig <- res[ res$padj < 0.1, ]
  goodRow <- resSig[complete.cases(resSig), ]
  geneSig <- rownames(goodRow)

  observedRegulation <- as.matrix(res$log2FoldChange)
  observedRegulation[observedRegulation > 0] <- 1
  observedRegulation[observedRegulation < 1] <- -1
  rownames(observedRegulation) <- rownames(res)

  if(dim(goodRow)[1] == 0) {
    resSig <- res[ res$pval < 0.1, ]
    goodRow <- resSig[complete.cases(resSig), ]
    geneSig <- rownames(goodRow)
  }

  return (list(cds=cds, res=res, geneSig=geneSig, resSig=resSig, observedRegulation=observedRegulation))
}

setwd(SAVE_DIR)
diffRes = differentialExpressedGene(countTable, condition)
write.csv(diffRes$geneSig, file = "diffRes.csv")
write.csv(diffRes$observedRegulation, file="observedRegulation.csv")
