#' Compute scaled or absolute population indices from a fitted model object, or from posterior samples.
#'
#' @param object A matrix or an object of class gam or brmsfit. If a matrix, columns should correspond to
#'               (posterior) samples of predictions and rows should match the rows in newdata.
#' @param timevar The name of the time variable in newdata over which an index should be computed.
#' @param byvar Name of grouping variable in newdata. The default is NULL in which case a single index is computed.
#'              If not null an index is computed for each unique value of the grouping variable.
#' @param newdata A data frame containing a time variable, any grouping variable, and any variables needed
#'                for predicting from object. The data frame will be supplied to the predict method for gam and brms models.
#'                If object is of class gam, any variables needed for prediction that are not available in
#'                newdata will be treated as constant when computing index.
#'                For objects of class brmsfit, all variables needed for predicting need to be available in newdata.
#'                If object is a matrix, the rows of newdata should correspond to the rows of the matrix.
#' @param type    Type of index to compute, one of "group", "global", "delta", "raw". If "group", indices for each group
#'                 are computed relative to the within group baseline. If "global", relative indices for each group are
#'                 computed relative to the global baseline. If "delta", indices for each group are computed relative to
#'                 the previous time point. If "raw", absolute (as opposed to relative) indices are computed.
#' @param weights  Weights for prediction points.
#' @param bweights Weights for the baseline, if different from target weights. If the argument is NULL and
#'                 there are multiple time points in the baseline, the average is taken (assumes newdata is balanced).
#' @param baseline A set of time points that should be used as baseline for indices of type "group" or "global".
#'                 The mean of the index over these time points will be one (see Knape 2023).
#'                 If missing, the first time point will be used as the baseline.
#' @param alpha A vector of alpha levels for computing confidence intervals.
#' @param nsamp Number of simulation samples to draw from gam objects. Defaults to NULL,
#'              in which case 1000 samples will be drawn.
#' @param ... Further arguments passed to predict functions.
#'
#'
#' @details
#' ## Warning
#' Large prediction tasks can require substantial memory. Use less rows in newdata,
#' and/or fewer simulation samples to reduce memory footprint.
#'
#' @return A data frame containing indices and their uncertainties.
#'
#' @examples
#' library(mgcv)
#' data(cuckoo)
#' \donttest{
#' # Simple model with abundance varying by year and latitude only.
#' gam_fit = gam(count ~ s(yr, lat), data = cuckoo, family = quasipoisson)
#'
#' # Compute index relative to first year at three example latitudes
#' nd = expand.grid(yr = unique(cuckoo$yr), lat = c(55,60,65))
#' index(gam_fit, time = "yr", newdata = nd)
#' }
#'
#' @export
#'
#' @references{Knape, J. (2023). Effects of choice of baseline on the uncertainty of population and biodiversity indices.
#'                               Environmental and Ecological Statistics, 30, 1--16. \doi{10.1007/s10651-022-00550-7}}
index = function(object, newdata, timevar, ...,  byvar = NULL,  type = "group", weights = NULL, bweights = NULL, baseline = NULL, alpha = c(.8, .95), nsamp = NULL) {
  type = match.arg(type, c("group", "global", "delta", "raw"))

  if (!is.character(timevar))
    stop(paste('Argument timevar should be a character string.'))

  if (!is.null(byvar)) {
    if (!is.character(byvar) | length(byvar) != 1)
      stop(paste('byvar argument should be a character string.'))
  }

  if (is.null(baseline) & (identical(type, "group") | identical(type, "global"))) {
    if (is.factor(newdata[[timevar]])) {
      baseline = levels(newdata[[timevar]])[1]
    } else {
      baseline = min(newdata[[timevar]])
    }
    if (!all(baseline %in% newdata[[timevar]]))
      stop("Baseline must contain time points that has observations.")
  }

  checkFamily(object)

  checkModel(object)

  newdata = checkNewdata(object, newdata, timevar)

  if (!(any(timevar == colnames(newdata)))) {
    stop(paste("Variable", timevar, "not found in newdata."))
  }
  if (!is.null(byvar)) {
    if (!(any(byvar == colnames(newdata)))) {
      stop(paste("Variable", byvar, "not found in newdata."))
    }
  }

  ##########################################################
  # Create index linking predictions to index and baselines
  ###########################################################
  if (!is.null(byvar)) {
    groupInd = intIndexDF(newdata[byvar])
  } else {
    groupInd = data.frame(key = rep(1, nrow(newdata)))
  }

  if (identical(type, "group")) {
    baseKey = groupInd[["key"]] * pmin(1, match(newdata[[timevar]], baseline)) # Integer from 1 to nbaselines, NA is skipped.
    baseInd = groupInd[["key"]]
  #  browser(), Fix if no observations matching baseline
    baseInd[!baseInd %in% baseKey] = NA
    nbase = length(baseline)
  }

  # Global indices
  if (identical(type, "global")) {
    baseKey = pmin(1, match(newdata[[timevar]], baseline))
    baseInd = rep(1, nrow(newdata))
    nbase = length(baseline)
  }

  # Delta indices
  if (identical(type, "delta")) {
    sortTimes = sort(unique(newdata[[timevar]]))
    if (length(sortTimes) < 2) {
      stop("Only one time point found, cannot estimate change.")
    }
    baseKey = (length(sortTimes) - 1)*(groupInd[["key"]] - 1) + match(newdata[[timevar]], sortTimes[-which.max(sortTimes)])
    baseInd = (length(sortTimes) - 1)*(groupInd[["key"]] - 1) + match(newdata[[timevar]], sortTimes[-1])
    nbase = 1
    if (!is.null(baseline)) {
      writeLines("baseline argument not used with type = \"delta\", ignored")
    }
  }

  if (identical(type, "raw")) {
    # TODO: Should check that all variables have been supplied to newdata, as constant effects don't cancel.
    baseInd = NULL
    baseKey = NULL
    if (!is.null(baseline)) {
      writeLines("baseline argument not used with type = \"raw\", ignored")
    }
    nbase = 0
  }

  if (identical(type, "custom")) {
    stop("Not implemented.")
    nbase="A"
  }

  ######################
  # Weights
  ######################

  if (is.null(weights)) {
    weights = rep(1, nrow(newdata))
  } else if (length(weights) == 1) {
    weights = rep(weights, nrow(newdata))
  } else if (length(weights) != nrow(newdata)) {
    stop("Length of weights does not match number of rows in newdata.")
  }
  if (!is.numeric(weights)) {
    stop("weights need to be numeric.")
  }

  # Baseline weights, default to 1/nbase * weights if not supplied.
  if (is.null(bweights)) {
    if (nbase > 1) {
      bweights = 1 / nbase * weights
    } else {
      bweights = weights
    }
  }

  if (length(bweights) == 1) {
    bweights = rep(bweights, nrow(newdata))
  } else if (length(bweights) != nrow(newdata)) {
    stop("Length of weights does not match number of rows in newdata.")
  }
  if (!is.numeric(bweights)) {
    stop("bweights need to be numeric.")
  }


  dataInd = newdata[, c(byvar, timevar), drop = FALSE]
  dataInd$baseInd = baseInd
  dataInd = intIndexDF(dataInd)
  dataInd$baseKey = baseKey
  attr(dataInd, "nbase") = nbase

  # Check for balance of weights
  if (length(unique(groupInd$key)) < 2) {
    tweights = tapply(weights, dataInd$key, sum)
    if (length(unique(tweights)) != 1) {
      warning("Total weights or number of prediction points differ across time points.")
    }
  } else {
    tweights = tapply(weights, list(dataInd$key, groupInd$key), sum)
    if (!all(sapply(apply(tweights, 2, unique, simplify = FALSE), length) == 2)) {
      warning("Total weights or number of prediction points differ across time points in at least one group.")
    }
  }


  opt = getSampOptions(object, nsamp, ...)

  extractor = getExtractor(object, opt)

  ind = computeind(extractor, dataInd, newdata = newdata, weights = weights, bweights = bweights, alpha = alpha, type = type, opt, ...)
  if (!is.null(baseline)) {
    attr(ind, "baseline") = baseline
  }
  attr(ind, "type") = type
  attr(ind, "nsamp") = nsamp
  #attr(ind, "constants") = constants
  attr(ind, "timevar") = timevar
  attr(ind, "byvar") = byvar
  attr(ind, "newdata") = newdata
  ind
}

checkFamily = function(object) {
  UseMethod("checkFamily")
}

checkModel = function(object) {
  UseMethod("checkModel")
}

checkNewdata = function(object, ...) {
  UseMethod("checkNewdata")
}

# If predictions can be constructed directly from object, then the method can simply return it, see matrix-method
getExtractor = function(object, ...)  {
  UseMethod("getExtractor")
}

getSampOptions = function(object, nsamp, ...) {
  UseMethod("getSampOptions")
}

# Note that repeat calls to getMu need to use the same random draws.
getMu = function(object, newdata, ...)  {
  UseMethod("getMu")
}



