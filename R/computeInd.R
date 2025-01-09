computeind = function(object, rowIndex, newdata, weights, bweights, alpha, type, opt, ...) {
  stopifnot(all(c(nrow(rowIndex), length(weights)) == nrow(newdata)))

  maxMatSize = 5e7L
  nsamp = opt$nsamp

  nd = nrow(newdata)
  size.mu = nd * nsamp
  nc = nsamp + identical(opt$pointestimate, "col1")
  exists.mu = FALSE

  ###################
  # Compute baseline
  ###################
  if (!identical(type, "raw")) {
    baseMat = matrix(NA, nrow = max(rowIndex[["baseKey"]], na.rm = TRUE), ncol = nc)
    maxBase = matrix(-Inf, nrow = nrow(baseMat), ncol = nc)
  }

  ###################
  # Compute index
  ###################

  outInd = rowIndex
  outInd[["baseKey"]] = NULL
  outInd = na.omit(unique(outInd))
  if (!identical(range(outInd[["key"]]), c(1L, nrow(outInd))))
    stop("id error")

  indMat = matrix(NA, nrow = nrow(outInd), ncol = nc)

  maxInd = matrix(-Inf, nrow = nrow(indMat), ncol = nc)

  batchIndMat = matrix(1:nrow(rowIndex))
  n.batch = ncol(batchIndMat)
  for (j in 1:n.batch) {
    batch.ind =  as.integer(na.omit(batchIndMat[,j]))
    muBatch = getMu(object, newdata[batch.ind, ], ...)
    stopifnot(ncol(muBatch) == nc)
    if (!is.null(rowIndex[["baseKey"]])) {
      maxBase = updColMaxBy(muBatch, maxBase, rowIndex[["baseKey"]][batch.ind], skip = is.na(rowIndex[["baseKey"]][batch.ind]))
      baseMat = updBaseMat(baseMat, muBatch, maxBase, rowIndex[["baseKey"]][batch.ind], weights = bweights[batch.ind], skip = is.na(rowIndex[["baseKey"]][batch.ind]))
    } else  {
      maxBase = baseMat = NULL
    }
    maxInd = updColMaxBy(muBatch, maxInd, rowIndex[["key"]][batch.ind], skip = is.na(rowIndex[["key"]][batch.ind]))
    indMat = updBaseMat(indMat, muBatch, maxInd, rowIndex[["key"]][batch.ind], weights = weights[batch.ind], skip = is.na(rowIndex[["key"]][batch.ind]))


  }
  out = indexStats(indMat, alpha = alpha, baseline = baseMat[outInd$baseInd, , drop = FALSE], pointestimate = opt$pointestimate)
  att.ci = attr(out, "ci")
  keep = setdiff(colnames(outInd), c("baseInd", "key"))
  out = as.data.frame(cbind(outInd[keep], out))
  attr(out, "ci") = att.ci
  out
}

