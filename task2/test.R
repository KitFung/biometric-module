
COUNT_PATH = 'raw_countTable.csv'
META_PATH = 'raw_meta.csv'

# processing raw data
countTable = read.csv(COUNT_PATH, header=TRUE, row.names = 1, check.names=FALSE)
metaMat = read.csv(META_PATH, header=TRUE, row.names = 1, check.names=FALSE)

bmi = lapply(as.vector(metaMat['体重指数', ]), as.character)
condition = res <- rep("condA", length(bmi))
condition[bmi <= 24] <- "condB"  # bmi <= 24


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
  resSig <- res[ res$padj < 0.1, ]
  goodRow <- resSig[complete.cases(resSig), ]
  geneSig <- rownames(goodRow)

  if(dim(goodRow)[1] == 0) {
    resSig <- res[ res$pval < 0.1, ]
    goodRow <- resSig[complete.cases(resSig), ]
    geneSig <- rownames(goodRow)
  }

  return (list(cds=cds, res=res, geneSig=geneSig, resSig=resSig))
}

#
diffRes = differentialExpressedGene(countTable, condition)
head( dres$res[ order(dres$res$pval), ] )
