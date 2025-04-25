#' @export
checkFamily.brmsfit = function(object) {
  family = family(object)
  if (!is.null(family$link_zi) | !is.null(family$link_hu)) {
    stop("Hurdle models and zero inflation not yet implemented/verified.")
  }
  if (!identical(family$link, "log")) {
    stop("Only log links allowed.")
  }
  family
}

#' @export
checkModel.brmsfit = function(object) {
}

#' @export
checkNewdata.brmsfit = function(object, newdata, ...) {
  newdata
}

#' @export
getExtractor.brmsfit = function(object, ...)  {
  object # If no batching, could instead supply newdata, and just return sample matrix.
}

#' @export
getSampOptions.brmsfit = function(object, nsamp, ...) {
  if (!is.null(nsamp)) {
    message("nsamp ignored for brms models, use argument ndraws to adjust number of samples")
  }
  dots = list(...)
  nsamp = brms::ndraws(object)
  if (!is.null(dots$ndraws)) {
    nsamp = dots$ndraws
  }
  if (!is.null(dots$draw_ids)) {
    nsamp = length(dots$draw_ids)
  }
  list(pointestimate = "mean", allow_batch = FALSE, nsamp = nsamp)
}

#' @export
getMu.brmsfit = function(object, newdata, ...) {
  # Prevent some options in ..., e.g. 'point_estimate', 'dpar', 'nlpar'?
  t(log(brms::posterior_epred(object, newdata = newdata, sort = FALSE, ...)))
}
