#' @export
checkFamily.gam = function(object) {
  family = object$family
  if (!identical(family$link, "log")) { # Should allow other links. Also allow transformed responses?
    stop("Only log links are currently allowed.")
  }
  family
}

#' @export
checkModel.gam = function(object) {
  if (!is.null(object$rank)) {
    if (object$rank < length(coef(object))) {
      warning('Model is rank deficient.')
    }
  }
}

#' @export
checkNewdata.gam = function(object, newdata, timevar, ...) {
  vars = all.vars(object$pred.formula)

  gvars = getGamCols(object)

  oix = sapply(gvars, \(x) {x[["type"]]}) == "offset"
  if (any(oix)) {
    offsvars = gvars[[which(oix)]][["vars"]]
  } else {
    offsvars = character(0)
  }
  constants = setdiff(vars, colnames(newdata))
  for (v in constants) {
    # If v an offset, warn.
    # First look in model frame
    nv = object$model[[v]][1]
    # If not there (because variable is transformed in formula), try var.summary
    if (is.null(nv)) {
      nv = object$var.summary[[v]][1] # MAY FAIL: var.summary may have levels that are deleted in model fit.
    }
    newdata[[v]] = nv
    if (v %in% offsvars) {
      warning(paste0("Missing offset ", v, " set to ", nv,", treated as constant for prediction."))
    } else {
      writeLines(paste0("Missing variable ", v, " set to ", nv,", treated as constant for prediction."))
    }
  }
  getInteractions = function(x) {
    if (any(timevar == x$vars) & length(x$vars) > 1)  {
      return(setdiff(x$vars, timevar))
    } else return(NULL)
  }
  ignoredInteractions = intersect(constants, unique(do.call(c, lapply(gvars, getInteractions))))
  if (length(ignoredInteractions) > 0) {
    warning(paste("Interactions between", timevar, "and", paste(ignoredInteractions, collapse = ", "), "are ignored.",
                  "Results may be misleading. Supply all variables included in interactions as newdata to change this behaviour."))
  }
  newdata
}

#' @export
getSampOptions.gam = function(object,nsamp, ...) {
  if (is.null(nsamp)) {
    nsamp = 1000
  }
  list(pointestimate = "col1", allow_batch = TRUE, nsamp = nsamp)
}

# Saves the vcov matrix so that it doesn't need to be recomputed every time predictions are needed.
#' @export
getExtractor.gam = function(object, opt, ...) {
  nsamp = opt$nsamp
  out = list(model = object)
  cf = coef(object)
  if (nsamp > 0) {
    vc  = mgcv::vcov.gam(object, unconditional = TRUE)
    #chvc = mgcv::mroot(vc, method = 'svd') # From help(predict.gam)
    chvc = mgcv::mroot(vc, method = 'chol') # Faster
    tol <- 100 * .Machine$double.eps
    chol.err <- max(abs(vc-chvc %*% t(chvc)))
    if (chol.err>tol) warning("Root of covariance matrix inaccurate, results may be unreliable.")
    cfn = cf + chvc %*% matrix(rnorm(ncol(chvc) * nsamp), nrow = ncol(chvc), ncol = nsamp)
    cf = cbind(cf, cfn)
  }
  out$cf = cf
  class(out) = "gam.extractor"
  out
}

#' @export
getMu.gam.extractor = function(object, newdata, ...) {
  X = predict(object$model, newdata = newdata, type = "lpmatrix", ...)
  mu = X %*% object$cf
  if (!identical(attr(X, "model.offset"), 0)) {
    mu = mu + rep(attr(X, "model.offset"), ncol(mu))
  }
  if (!identical(object$model$family$link, "log")) { # Not yet allowed by checkFamily...
    mu = log(object$family$linkinv(mu))
  }
  mu
}

#' @export
checkFamily.gamm = function(object) {
  checkFamily(object[["gam"]])
}

#' @export
checkModel.gamm = function(object) {
  # Check PQL etc?
  checkModel(object[["gam"]])
}

#' @export
checkNewdata.gamm = function(object, newdata, timevar, ...) {
  checkNewdata(object[["gam"]], newdata, timevar, ...)
}

#' @export
getSampOptions.gamm = function(object, nsamp, ...) {
  getSampOptions(object[["gam"]], nsamp, ...)
}

#' @export
getExtractor.gamm = function(object, opt, ...) {
  getExtractor(object[["gam"]], opt, ...)
}


# Should work with both gam and gam.prefit
getGamCols = function(object) {
  # Parametric terms
  terms = delete.response(object$pterms)
  labs = attr(terms, 'term.labels')
  assign = object$assign
  pvars = NULL
  if (any(assign == 0)) {
    stopifnot(attr(terms, 'intercept') == 1)
    pvars[["(Intercept)"]] = list(type = "intercept",
                                  vars = NULL,
                                  cols = which(assign == 0))
  }
  for (l in seq_len(length(labs))) {
    pvars[[labs[l]]] = list(type = "parametric", vars = all.vars(str2expression(labs[l])), cols = which(assign == l))
    # Simple check for mismatch, not bulletproof
    for (v in pvars[[labs[l]]]$vars) {
      if (any(!grepl(v, names(object$cmX)[pvars[[labs[l]]]$cols], fixed = TRUE)))
        stop("Column mismatch")
    }
  }
  # Smooth terms
  svars = vector(length(object$smooth), mode = "list")
  smoothTerms = lapply(object$smooth, `[`, c('term', 'label', 'by', 'first.para', 'last.para'))

  for (l in seq_len(length(svars))) {
    names(svars)[l] = smoothTerms[[l]]$label
    vars = all.vars(str2expression(smoothTerms[[l]]$term))
    if (!identical(smoothTerms[[l]]$by, "NA"))
      vars = c(vars, all.vars(str2expression(smoothTerms[[l]]$by)))
    svars[[l]] = list(type = "smooth",
                      vars = vars,
                      cols = smoothTerms[[l]]$first.para:smoothTerms[[l]]$last.para)
  }
  off = attr(terms, "offset")
  if (!is.null(off)) {
    ovars = list(list(type = "offset", vars = all.vars(attr(terms,"variables")[[off + 1]]), cols = NULL))
    names(ovars) = deparse(attr(terms,"variables")[[off + 1]])
  } else {
    ovars = NULL
  }
  vars = c(pvars, svars, ovars)

  if (!identical(as.integer(do.call(c, lapply(vars, `[[`, 'cols'))), seq_len(length(object$cmX))))
    stop("Missing columns")
  if (!setequal(unique(do.call(c, lapply(vars, `[[`, 'vars'))), all.vars(object$pred.formula)))
    stop("Variable mismatch")
  vars
}

