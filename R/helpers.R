
intIndexDF = function(df, return.all = TRUE) {
  stopifnot(is.data.frame(df))
  fac = do.call(interaction, c(df, list(drop = TRUE)))
  df$key = as.integer(factor(fac, levels = unique(fac)))
  if (!return.all)
    df = df["key"]
  df
}


indexStats = function(ind, alpha, log = FALSE, prefix = "", pointestimate = "col1", baseline = NULL) { # Assumes that first column is point estimate
  nr = nrow(ind)
  pointestimate = match.arg(pointestimate, c("col1", "mean", "median", "log_mean"))
  if (!is.null(baseline)) {
    ind = ind - baseline
    if (identical(pointestimate, "col1")) {
      bs = baseline[, 1]
      baseline = baseline[, -1, drop = FALSE]
    } else {
      bs = switch(pointestimate,
                col1 = bs,
                mean = log(rowMeans(exp(baseline))),
                median = apply(baseline, 1, median),
                log_mean = rowMeans(baseline))
    }
  } else{
    bs= NULL
  }
  if (identical(pointestimate, "col1")) {
    index = ind[,1]
    ind = ind[, -1, drop = FALSE]
  } else {
    index = switch(pointestimate,
                 col1 = index,
                 mean = log(rowMeans(exp(ind))),
                 median = apply(ind, 1, median),
                 log_mean = rowMeans(ind))
  }
  ind = na.omit(ind)

  probs = sort(c((1 - alpha)/2, 1 - (1 - alpha)/2))
  if (!is.null(alpha)) {
    if (nrow(ind > 0)) {
      ci = t(apply(ind, 1, stats::quantile, probs, names = FALSE))
    } else {
      ci = matrix(nrow = 0, ncol = length(probs))
    }
    levnames = formatC(100*probs, digits = 3, width = 1, format = "fg")
    colnames(ci) = paste0(prefix, 'ci_', levnames)
  } else {
    ci = matrix(nrow = nrow(ind), ncol = 0)
    levnames = character(0)
  }
  se_log = apply(ind, 1, sd)
  se = apply(exp(ind), 1, sd)
  ses = cbind(se = se, se_log = se_log)
  colnames(ses) = paste0(prefix, colnames(ses))
  out = matrix(nrow = nr, ncol = 1 + ncol(ci) + ncol(ses))
  if (log) {
    out[,1] = index
  } else {
    out[,1] = exp(index)
    ci = exp(ci)
  }
  if (!is.null(attr(ind, 'na.action')))
    out[-attr(ind, 'na.action'), 2:ncol(out)] = cbind(ci, ses)
  else
    out[, 2:ncol(out)] = cbind(ci, ses)
  colnames(out) = c("index", colnames(ci), colnames(ses))
  if (!is.null(bs)) {
    if (log) {
      out = cbind(out, baseline = bs)
    } else {
      out = cbind(out, baseline= exp(bs))
    }
  }

  if (!is.null(alpha)) {
    attr(out, 'ci') = data.frame(quantile = probs, name = levnames, id = c(1:length(alpha), length(alpha):1), lu = rep(c('l', 'u'), each = length(alpha)))
  }
  out
}


