# metaDat: matrix, row = sample, col = attribute
# usedAttribute = array of used attribute
# combinations: matrix, 1 row = 1 combination. Each columns represent the corresponded value in usedAttribute
# return
separateMetaGroup <- function(metaMat, usedAttribute, combinMat) {
  if(length(usedAttribute) != dim(combinMat)[2]) {
    stop("Length of usedAttribute not equal to the width of combinMat")
  }

  sampleCnt <- dim(metaMat)[1]
  combinCnt <- dim(combinMat)[1]
  usedAttributeCnt <- length(usedAttribute)
  res <- rep(1, sampleCnt)

  for (i in 1:combinCnt) {
    cond <- metaMat[, usedAttribute] == combinMat[i,]
    matchIdx <- cond
    if(class(cond) == "matrix")
      matchIdx <- rowSums(cond) == usedAttributeCnt
    res[matchIdx] <- 0;
  }

  return (factor(res))
}
