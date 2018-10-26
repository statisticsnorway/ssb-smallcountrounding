

# stackoverflow questions 30357330
pkgEnvSmallCountData <- new.env(parent=emptyenv())
#' Function that returns a dataset 
#'
#' @encoding UTF8
#'
#' @param dataset Name of data set within the SmallCountRounding  package  
#' @param path When non-NULL the data set is read from "path/dataset.RData"
#'
#' @return The dataset
#' @export
#' @importFrom utils data
#' 
#' @note Currently, the function returns the same datasets as is included in the package easySdcTable.
#'
#' @examples
#'  z  <- SmallCountData("sosialFiktiv")
#'
SmallCountData <- function(dataset, path = NULL) {
  if (!is.null(path)) {
    filename <- paste(path, dataset, ".RData", sep = "")
    return(get(load(filename, envir = environment())))
  }
  if (!exists(dataset, pkgEnvSmallCountData))
    data(list = dataset, package = "SmallCountRounding", envir = pkgEnvSmallCountData)
  return(pkgEnvSmallCountData[[dataset]])
  return(NULL)
}

#' Fictitious datasets used in the examples.
#' 
#' The most comprehensive dataset, \code{sosialFiktiv}, contains three dimensions. The first dimension is 'region' which is grouped in two ways, 'fylke' and  'kostragr'. The other two are 'hovedint' and 'mnd'. In 'mnd2' two of the three categories in 'mnd' are merged.
#' The other datasets (\code{z1}, \code{z1w}, \code{z2}, \code{z2w}, \code{z3}, \code{z3w}, \code{z3wb}) are smaller subdatasets.
#' Datasets marked with '\code{w}' are unstacked and several variables are holding counts.
#'
#' @docType data
#' @keywords datasets internal
#' @name sosialFiktiv
NULL

#' @rdname sosialFiktiv
#' @name z1
NULL

#' @rdname sosialFiktiv
#' @name z1micro
NULL

#' @rdname sosialFiktiv
#' @name z1w
NULL

#' @rdname sosialFiktiv
#' @name z2
NULL

#' @rdname sosialFiktiv
#' @name z2w
NULL

#' @rdname sosialFiktiv
#' @name z3
NULL

#' @rdname sosialFiktiv
#' @name z3w
NULL

#' @rdname sosialFiktiv
#' @name z3wb
NULL