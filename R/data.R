#' Counts of common cuckoo from the Swedish Bird Survey
#'
#' Counts per survey route and year between 2000 and 2021 of common cuckoos from the
#' Swedish Bird Survey.
#'
#' @format ## `cuckoo`
#' A data frame with 6 columns.
#' \describe{
#'   \item{count}{Number of individuals counted.}
#'   \item{yr}{Year}
#'   \item{route}{Id of the survey route.}
#'   \item{county}{County where route is located.}
#'   \item{lon}{Longitude}
#'   \item{lon}{Latitude}
#' }
#' @source \doi{10.15468/hd6w0r}
"cuckoo"


#' Map of Swedish counties
#'
#' Polygons of Swedish counties and their areas.
#'
#' @format ## `cuckoo`
#' An sf data frame with 3 columns.
#' \describe{
#'   \item{county}{Name of county.}
#'   \item{area}{Area of county in km^2.}
#'   \item{geometry}{sf geometry column}
#' }
#' @source <https://www.scb.se/hitta-statistik/regional-statistik-och-kartor/regionala-indelningar/digitala-granser/>
"swe_map"
