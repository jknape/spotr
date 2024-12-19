#' @export
checkFamily.matrix = function(object) {
  if (!is.numeric(object))
    stop("Samples should be numeric.")
  if (any(object < 0, na.rm = TRUE)) {
    stop("Samples need to be positive.")
  }
}

#' @export
checkModel.matrix = function(object) {
  stopifnot(is.matrix(object))
}

#' @export
checkNewdata.matrix = function(object, newdata, ...) {
  if (!identical(nrow(object), nrow(newdata)))
    stop("Number of rows in sample matrix need to match number of rows in newdata.")
  newdata
}

#' @export
getExtractor.matrix = function(object, ...)  {
  log(object)
}

#' @export
getSampOptions.matrix = function(object, nsamp, ...) {
  if (!is.null(nsamp)) {
    writeLines("Using all columns in matrix, nsamp ignored.")
  }
  list(pointestimate = "mean", allow_batch = FALSE, nsamp = ncol(object))
}

#' @export
getMu.matrix = function(object, ...) {
  object
}

