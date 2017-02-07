
COUNT_PATH = 'raw_countTable.csv'
META_PATH = 'raw_meta.csv'

# processing raw data
metaMat = read.csv(META_PATH, header=TRUE, row.names = 1, check.names=FALSE)
countTable = read.csv(COUNT_PATH, header=TRUE, row.names = 1, check.names=FALSE)
countTable = countTable[rowSums(countTable) != 0, ]

bmi = lapply(as.vector(metaMat['体重指数', ]), as.character)
condition = rep("condB", length(bmi))
condition[bmi <= 24] <- "condA"  # bmi <= 24


# ##IF the meta data is in format columns is attribute, row is sample
# ## and the condition is just equal to something, can use this.
# condition <- separateMetaGroup(metaMat, usedAttribute, combinMat)

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

diffRes = differentialExpressedGene(countTable, condition)
write.csv(diffRes$geneSig, file = "diffRes.csv")
write.csv(diffRes$observedRegulation, file="observedRegulation.csv")
countTable
condition <- as.matrix(condition)
rownames(condition) <- colnames(metaMat)
write.csv(condition, file = "condition.csv")
