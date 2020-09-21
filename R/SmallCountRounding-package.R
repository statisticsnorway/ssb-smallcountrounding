#' Small Count Rounding of Tabular Data 
#' 
#' A statistical disclosure control tool to protect frequency tables in cases where small values are sensitive.
#' The main function, \code{\link{PLSrounding}}, performs small count rounding of necessary inner cells (Heldal, 2017) 
#' so that all small frequencies of cross-classifications to be published (publishable cells) are rounded. 
#' This is equivalent to changing micro data since frequencies of unique combinations are changed. 
#' Thus, additivity and consistency are guaranteed.
#' This is performed by an algorithm inspired by partial least squares regression (Langsrud and Heldal, 2018). \cr
#' 
#' @name SmallCountRounding-package
#' @aliases SmallCountRounding
#' @docType package
#' 
#' @references 
#' Heldal, J. (2017): \dQuote{The European Census Hub 2011 Hypercubes - Norwegian SDC Experiences}. 
#' In: \emph{Work Session on Statistical Data Confidentiality}, 
#' Skopje, The former Yugoslav Republic of Macedonia, September 20-22 , 2017.
#' 
#' Langsrud, Ø. and Heldal, J. (2018): \dQuote{An Algorithm for Small Count Rounding of Tabular Data}. 
#' Presented at: \emph{Privacy in statistical databases}, Valencia, Spain. September 26-28, 2018.
#' \url{https://www.researchgate.net/publication/327768398_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data}
#' 
#' @encoding UTF8
NULL