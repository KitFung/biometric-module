args = commandArgs(TRUE)
TRAIN_DATA_PATH = args[1]
TRAIN_LABEL_PATH = args[2]
TEST_DATA_PATH = args[3]
SAVE_DIR = file.path(args[4], "classify")


library(YuGene)
library(pamr)
library(siggenes)
library(ROCR)

findDiffGenes <- function(sdat.f, clus.f, unique_label, pvalth=0.01) {
  diffGenes <- NULL
  for(cl in unique_label) {
    temp <- clus.f
    temp[clus.f==cl] <- 0
    temp[clus.f!=cl] <- 1
    res.sam <- sam(sdat.f, c(t(temp)), gene.names=rownames(sdat.f), rand=123456)
    temp.diffGenes <- names(res.sam@q.value)[res.sam@q.value <= pvalth]
    diffGenes <- c(diffGenes, temp.diffGenes)
  }
  diffGenes <- unique(diffGenes)
  return(diffGenes)
}

filterDiffGenes <- function(sdat.f, clus.f, diffGenes, unique_label, aucth=0.9) {
  sigdat <- sdat.f[diffGenes, ]
  sigGeneIds <- rownames(sigdat)
  compareTwoClus <- function(sigdat, clus.f, clu=1) {

    temp <- clus.f
    temp[clus.f==clu] <- 0
    temp[clus.f!=clu] <- 1

    temp2 <- temp
    temp2[temp==0] <- 1
    temp2[temp==1] <- 0
    sigdat <- sigdat[, rownames(temp)]

    preds1 <- apply(sigdat, 1, prediction, labels = temp)
    preds2 <- apply(sigdat, 1, prediction, labels = temp2)

    aucs1 <- unlist(lapply(preds1, function(x) {
      (performance(x, measure = "auc"))@y.values[[1]]
    }))
    aucs2 <- unlist(lapply(preds2, function(x) {
      (performance(x, measure = "auc"))@y.values[[1]]
    }))
    aucs <- rowMax(cbind(aucs1, aucs2))

    diffGenes.f <- rownames(sigdat)[which(aucs >= aucth)]
    return(diffGenes.f)
  }
  df <- NULL
  for(i in unique_label) {
    df <- c(df, compareTwoClus(sigdat, clus.f,i))
  }
  diffGenes.f <- unique(df)
  return(diffGenes.f)
}


getCentroids <- function (fit, data, threshold) {

    genenames <- data$genenames[fit$gene.subset]
    x <- data$x[fit$gene.subset, fit$sample.subset]
    clabs <- colnames(fit$centroids)
    scen <- pamr.predict(fit, data$x, threshold = threshold,
        type = "cent")
    dif <- (scen - fit$centroid.overall)/fit$sd
    if (!is.null(fit$y)) {
        nc <- length(unique(fit$y))
    }
    if (is.null(fit$y)) {
        nc <- ncol(fit$proby)
    }
    o <- drop(abs(dif) %*% rep(1, nc)) > 0
    d <- dif[o, ]
    nd <- sum(o)
    genenames <- genenames[o]
    xx <- x[o, ]
    oo <- order(apply(abs(d), 1, max))
    d <- d[oo, ]
    genenames <- genenames[oo]
    return(d)
}

buildClassifier <- function(sigMat, clus.f, unique_label, nfold=10, nboot=100) {
    dat <- list()
    dat$x <- sigMat[, rownames(clus.f)]
    for(i in unique_label) {
      dat$y[clus.f==i] <- i
    }

    dat$y <- factor(dat$y)
    dat$geneid <- rownames(dat$x)
    dat$genenames <- rownames(dat$x)
    pam.rslt <- pamr.train(data=dat)
    pam.cv.rslt.l <- list()
    for(i in 1:nboot) {
        pam.cv.rslt <- pamr.cv(fit=pam.rslt, data=dat, nfold=nfold)
        pam.cv.rslt.l[[i]] <- pam.cv.rslt
    }
    err <- t(sapply(1:length(pam.cv.rslt.l), function(x)
        pam.cv.rslt.l[[x]]$error, simplify=TRUE))
    ngenes <- c(sapply(1:(length(pam.rslt$threshold)-1), function(x)
        nrow(pamr.listgenes(pam.rslt, dat, pam.rslt$threshold[x]))
    ), 0)
    colnames(err) <- ngenes
    thresh <- pam.rslt$threshold[18]

    signature <- (pamr.listgenes(pam.rslt, dat, thresh))[, "id"]

    cents <- getCentroids(pam.rslt, dat, thresh)
    cents <- cents[signature, ]

    return(list(signature=signature, pam.rslt=pam.rslt, thresh=thresh,
        err=err, cents=cents))
}

pamClassify <- function(datsel, signature, pam.rslt, thresh, unique_label,　postRth=1) {
    pred <- pamr.predict(pam.rslt, datsel, thresh, type="posterior")
    maxr <- apply(pred, 1, max)
    postR <- maxr/(1-maxr)

    sel.samples.f <- names(postR)[which(postR >= postRth)]

    clu.pred <- apply(pred, 1, which.max)
    clu.pred <- sort(clu.pred)
    clu.pred.reord<-NULL
    for(cl in unique_label) {
        temp <- names(sort(postR[names(clu.pred[clu.pred==cl])]))
        clu.pred.reord <- c(clu.pred.reord, temp)
    }
    clu.pred <- clu.pred[clu.pred.reord]

    nam.ord <- names(clu.pred)

    sdat.sig <- datsel[signature, nam.ord]
    gclu.f <- hclust(as.dist(1-cor(t(sdat.sig))), method="complete")

    return(list(sdat.sig=sdat.sig, pred=pred, clu.pred=clu.pred,
        nam.ord=nam.ord, gclu.f=gclu.f))
}


dir.create(SAVE_DIR, showWarnings = FALSE)
do.call(file.remove,list(list.files(path=SAVE_DIR, pattern="*", full.names=T)))

# load data
train_data <- as.matrix(
            read.table(TRAIN_DATA_PATH, header=TRUE, sep = "\t",
                row.names = 1,
                check.names=FALSE,
                as.is=TRUE))
storage.mode(train_data) <- "numeric"
train_label <- read.csv(TRAIN_LABEL_PATH, header=TRUE, row.names=1)

test_data <- as.matrix(
            read.table(TEST_DATA_PATH, header=TRUE, sep = "\t",
                row.names = 1,
                check.names=FALSE,
                as.is=TRUE))
storage.mode(test_data) <- "numeric"

unique_label <- unique(unlist(c(train_label)))

# normalize the train and test data
no_zero_train_data <- train_data[which(rowSums(train_data==0) != ncol(train_data)), ]
train_mat <- YuGene(no_zero_train_data)
test_mat <- YuGene(test_data)

# filtering sample
nofilterdata <- intersect(rownames(train_mat), rownames(test_mat))

dfg1 <- findDiffGenes(train_mat, train_label, unique_label)
diffGenes.f <- intersect(dfg1, rownames(test_mat))
if(length(diffGenes.f) < 5)
  diffGenes.f <- nofilterdata

diffGenes.f <- filterDiffGenes(train_mat, train_label, diffGenes.f, unique_label)
if(length(diffGenes.f) < 5)
  diffGenes.f <- nofilterdata

sigMat <- train_mat[diffGenes.f, rownames(train_label)]

# build classifier
classifier <- buildClassifier(sigMat, train_label, unique_label)

signature <- classifier$signature
pam.rslt <- classifier$pam.rslt
thresh <- classifier$thresh
err <- classifier$err

# do classify
datsel <- test_mat[diffGenes.f, rownames(train_label)]
pamcl <- pamClassify(datsel, signature, pam.rslt, thresh, unique_label,　postRth=1)
sdat.sig <- pamcl$sdat.sig
pred <- pamcl$pred
clu.pred <- pamcl$clu.pred
nam.ord <- pamcl$nam.ord
gclu.f <- pamcl$gclu.f

# save result
setwd(SAVE_DIR)
write.csv(clu.pred, file="classifed_result.csv")
write.csv(diffGenes.f, file="diffGenes.f.csv")
write.csv(thresh, file="thresh.csv")
write.csv(pam.rslt$centroids, file="centroids.csv")
write.csv(pam.rslt$centroid.overall, file="centroid.overall.csv")
