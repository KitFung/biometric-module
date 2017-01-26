args = commandArgs(TRUE)
DATA_PATH = args[1]
SAVE_DIR = file.path(args[2], "cluster")
MAX_K = 10

library(YuGene)
library(Biobase)
library(cluster)
library(siggenes)
library(ConsensusClusterPlus)
library(AnnotationDbi)
library(pamr)
library(gplots)
library(ROCR)

source("figures.R")

compGapStats <- function(ge.CRC, ntops=c(1, 2, 3)*100, K.max=6, nboot=100) {
    MAD <- apply(ge.CRC, 1, mad)
    ords <- names(sort(MAD, decreasing=TRUE))

    if(length(ords) < 300) {
      ntops <- c(2, 3, length(ords));
    }

    fun <- function(x, k) {
        y <- hclust(as.dist(1-cor(t(x))), method="average")
        clus <- cutree(y, k=k);
        return(list(cluster=clus))
    }
    gaps <- list()
    for(i in 1:length(ntops)) {
        ntop <- ntops[i]
        sdat <- ge.CRC[ords[1:ntop], ]
        sdat <- sweep(sdat,1, apply(sdat,1,median))
        gaps[[i]] <- clusGap(t(sdat), FUNcluster=fun, K.max=K.max, B=nboot)
    }
    gapsmat <- matrix(0, K.max, length(ntops))
    gapsSE <- matrix(0, K.max, length(ntops))
    for(i in 1:length(ntops)) {
        gapsmat[, i] <- gaps[[i]]$Tab[, 3]
        gapsSE[, i] <- gaps[[i]]$Tab[, 4]
    }

    colnames(gapsmat) <- colnames(gapsSE) <- ntops
    rownames(gapsmat) <- rownames(gapsSE) <- 1:K.max
    return(list(gapsmat=gapsmat, gapsSE=gapsSE))
}


selTopVarGenes <- function(ge.CRC, MADth=0.5) {
    MAD <- apply(ge.CRC, 1, mad)
    sdat <- ge.CRC[which(MAD > MADth), ]
    sdat = sweep(sdat,1, apply(sdat,1,median))
    return(sdat)
}

conClust <- function(sdat, bestK=3, maxK=12, reps=1000, savepath=".") {
    res <- ConsensusClusterPlus(sdat, maxK=maxK, reps=reps, pItem=0.98,
        innerLinkage="average", finalLinkage="average", pFeature=1,
        title=savepath, clusterAlg="hc", distance="pearson", plot="pdf")
    areaK <- res[["areaK"]]
    res[["areaK"]] <- NULL
    icl <- calcICL(res, title=savepath, plot="pdf")

    con <- icl[[2]][icl[[2]][, "k"]==bestK, ]
    conMat <- NULL
    for(i in 1:bestK) {
        conMat <- cbind(conMat, con[con[, 2]==i, 4])
    }
    rownames(conMat) <- con[con[, 2]==1, 3]
    clus <- apply(conMat, 1, which.max)
    return(clus)
}

filterSamples <- function(sdat, clus) {
    sdat.ord <- sdat[, names(clus)]
    silh <- silhouette(clus, as.dist(1-cor(sdat.ord)))
    excl <- names(clus)[which(silh[, "sil_width"] < 0)]
    sdat.f <- sdat.ord[, setdiff(colnames(sdat.ord), excl)]
    clus.f <- clus[which(names(clus) %in% colnames(sdat.f))]
    clus.f <- sort(clus.f)

    return(list(sdat.f=sdat.f, clus.f=clus.f, silh=silh))
}


dir.create(SAVE_DIR, showWarnings = FALSE)
do.call(file.remove,list(list.files(path=SAVE_DIR, pattern="*", full.names=T)))

exprs <- as.matrix(
            read.table(DATA_PATH, header=TRUE, sep = "\t",
                row.names = 1,
                check.names=FALSE,
                as.is=TRUE))
storage.mode(exprs) <- "numeric"

# Preprocessing
no_zero_exprs <- exprs[which(rowSums(exprs==0) != ncol(exprs)), ]
mat <- YuGene(no_zero_exprs)

# Compute Gap statistic
gaps <- compGapStats(mat, nboot=10)
gapsmat <- gaps$gapsmat # gap value
gapsSE <- gaps$gapsSE # standard error of gap

# use the best K in top 300
bestK <-  which(max(gapsmat[, 3]) == gapsmat[, 3])

sdat <- selTopVarGenes(mat, MADth=0.01) # 0.5

clus <- conClust(sdat, bestK=bestK, savepath=SAVE_DIR)

samp.f <- filterSamples(sdat, clus)
sdat.f <- samp.f$sdat.f
clus.f <- samp.f$clus.f
silh <- samp.f$silh

cat("Done process\n")

setwd(SAVE_DIR)
write.csv(clus, file="consesClusterResult.csv")
##
pdf(file="gap.pdf", width=8, height=8)
figGAP(gapsmat, gapsSE)
graphics.off()
##
pdf(file="silh.pdf", width=8, height=8)
figSilh(silh)
graphics.off()
##

cat("Done all\n")
