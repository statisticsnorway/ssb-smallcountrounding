#' PLS inspired rounding
#' 
#'  Small count rounding of necessary inner cells are performed so that all small frequencies of cross-classifications to be published 
#' (publishable cells) are rounded. The publishable cells can be defined from a model formula, hierarchies or automatically from data.
#' 
#' This function is a user-friendly wrapper for \code{RoundViaDummy} with data frame output and with computed summary of the results.
#' See \code{\link{RoundViaDummy}} for more details.
#'
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (inner cells frequencies)
#' @param roundBase Rounding base
#' @param hierarchies List of hierarchies
#' @param formula Model formula defining publishable cells
#' @param ... Further parameters sent to \code{RoundViaDummy}  
#'
#' @return Output is a four-element list with class attribute "PLSrounded" (to ensure informative printing).
#'    \item{inner}{Data frame corresponding to input data with the main dimensional variables and with cell 
#'                frequencies (original, rounded, difference).}
#'    \item{publish}{Data frame of publishable data with the main dimensional variables and with cell frequencies 
#'                   (original, rounded, difference).}
#'    \item{metrics}{A named character vector of various statistics calculated from the two output data frames 
#'    ("\code{inner_}" used to distinguish). See examples below and the function \code{\link{HDutility}}.}
#'    \item{freqTable}{Matrix of frequencies of cell frequencies and absolute differences.
#'    For example, row "\code{rounded}" and column "\code{pub.4+}" is the number of rounded 
#'    inner cell frequencies greater than or equal to \code{4}.}
#' @export
#'
#' @examples
#' # Small example data set
#' z <- SmallCountData("e6")
#' print(z)
#' 
#' # Publishable cells by formula interface
#' a <- PLSrounding(z, "freq", roundBase = 5,  formula = ~geo + eu + year)
#' print(a)
#' print(a$inner)
#' print(a$publish)
#' print(a$metrics)
#' print(a$freqTable)
#' 
#' # Recalculation of maxdiff, HDutility, meanAbsDiff and rootMeanSquare
#' max(abs(a$publish[, "difference"]))
#' HDutility(a$publish[, "original"], a$publish[, "rounded"])
#' mean(abs(a$publish[, "difference"]))
#' sqrt(mean((a$publish[, "difference"])^2))
#' 
#' # Four lines below produce equivalent results 
#' # Ordering of rows can be different
#' PLSrounding(z, "freq")
#' PLSrounding(z, "freq", formula = ~eu * year + geo * year)
#' PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eHrc"))
#' PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"))
#' 
#' # Define publishable cells differently by making use of formula interface
#' PLSrounding(z, "freq", formula = ~eu * year + geo)
#' 
#' # Define publishable cells differently by making use of hierarchy interface
#' eHrc2 <- list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019"))
#' PLSrounding(z, "freq", hierarchies = eHrc2)
#' 
#' # Package sdcHierarchies can be used to create hierarchies. 
#' # The small example code below works if this package is available. 
#' if (require(sdcHierarchies)) {
#'   z2 <- cbind(geo = c("11", "21", "22"), z[, 3:4], stringsAsFactors = FALSE)
#'   h2 <- list(
#'     geo = hier_compute(inp = unique(z2$geo), dim_spec = c(1, 1), root = "Tot", as = "df"),
#'     year = hier_convert(hier_create(root = "Total", nodes = c("2018", "2019")), as = "df"))
#'   PLSrounding(z2, "freq", hierarchies = h2)
#' }
PLSrounding <- function(data, freqVar, roundBase = 3, hierarchies = NULL, formula = NULL, ...) {
  z <- RoundViaDummy(data = data, freqVar = freqVar, formula = formula, roundBase = roundBase, hierarchies = hierarchies, ...)
  
  
  MakeDifference <- function(x) cbind(x, difference = x[, "rounded"] - x[, "original"])
  
  z$yInner <- MakeDifference(z$yInner)
  z$yPublish <- MakeDifference(z$yPublish)
  
  maxdiff <- max(abs(z$yPublish[, "difference"]))
  
  inner_HDutility <- HDutility(z$yInner[, "original"], z$yInner[, "rounded"])
  publish_HDutility <- HDutility(z$yPublish[, "original"], z$yPublish[, "rounded"])
  
  inner_meanAbsDiff <- mean(abs(z$yInner[, "difference"]))
  publish_meanAbsDiff <- mean(abs(z$yPublish[, "difference"]))
  
  inner_rootMeanSquare <- sqrt(mean((z$yInner[, "difference"])^2))
  publish_rootMeanSquare <- sqrt(mean((z$yPublish[, "difference"])^2))
  
  freqTable <- cbind(TabCat(z$yInner, roundBase, "inn."), TabCat(z$yPublish, roundBase, "pub."))
  
  freqVarName <- names(data[1, freqVar, drop = FALSE])
  
  out <- NULL
  
  if (!is.null(z$crossTable)) {
    cNames <- colnames(data)[colnames(data) %in% colnames(z$crossTable)]
    out$inner <- cbind(data[, cNames, drop = FALSE], z$yInner)
    out$publish <- cbind(as.data.frame(z$crossTable[, cNames, drop = FALSE]), z$yPublish)
    rownames(out$publish) <- NULL
  } else {
    out$inner <- as.data.frame(z$yInner)
    out$publish <- as.data.frame(z$yPublish)
  }
  
  if (any(duplicated(colnames(out$inner)))) {
    warning(paste("Duplicated colnames in output:", paste(colnames(out$inner)[duplicated(colnames(out$inner))], collapse = ", ")))
  }
  
  
  out$metrics <- c(roundBase = roundBase, maxdiff = maxdiff, inner_HDutility = inner_HDutility, HDutility = publish_HDutility, 
                   inner_meanAbsDiff = inner_meanAbsDiff, meanAbsDiff = publish_meanAbsDiff, 
                   inner_rootMeanSquare = inner_rootMeanSquare, rootMeanSquare = publish_rootMeanSquare)
  out$freqTable <- freqTable
  out
  return(structure(out, class = "PLSrounded"))
}



#' Print method for PLSrounded
#'
#' @param x PLSrounded object 
#' @param digits positive integer.  Minimum number of significant digits to be used for printing most numbers.
#' @param \dots further arguments sent to the underlying
#'
#' @return Invisibly returns the original object.
#' @keywords print
#' @export
print.PLSrounded <- function(x, digits = max(getOption("digits") - 3, 3), ...) {
  b <- x$metrics["roundBase"]
  cat("\nPLSrounding summary:  \n\n")
  metricsToPrint <- x$metrics[c("maxdiff", "HDutility", "meanAbsDiff", "rootMeanSquare")]
  print(factor(round(metricsToPrint, 4)), max.levels = 0, digits = digits, ...)
  cat("\nFrequencies of cell frequencies and absolute differences:  \n\n")
  print.table(x$freqTable, zero.print = ".", digits = digits, ...)
  cat("\n")
  invisible(x)
}


TabCat <- function(z, b, s) {
  x <- rbind(original = table(ToCat(z[, "original"], b)), 
             rounded = table(ToCat(z[, "rounded"], b)), 
             absDiff = table(ToCat(abs(z[, "difference"]), b)))
  x <- x[, colnames(x) != "0.5", drop = FALSE]
  x <- cbind(x, all = NROW(z))
  colnames(x) <- paste(s, colnames(x), sep = "")
  x
}


ToCat <- function(x, b) {
  x[x > b] <- b + 1
  x[x > 0 & x < b] <- 1 - as.numeric(b == 1)/2
  x <- factor(x, levels = c(0, 1 - as.numeric(b == 1)/2, b, b + 1))
  levels(x)[4] <- paste(levels(x)[4], "+", sep = "")
  if (b > 2) 
    levels(x)[2] <- paste(1, as.character(b - 1), sep = "-")
  x
}


#' @rdname HDutility 
#' @export
HD <- function(f, g){ 
  sqrt(sum((sqrt(f) - sqrt(g))^2)/2)
}


#' Hellinger Distance (Utility)
#' 
#' Hellinger distance (\code{HD}) and a related utility measure (\code{HDutility})
#' described in the reference below.
#' The utility measure is made to be bounded between 0 and 1.
#' 
#' HD is defined as "\code{sqrt(sum((sqrt(f) - sqrt(g))^2)/2)}" and 
#' HDutility  is defined as "\code{1 - HD(f, g)/sqrt(sum(f))}".
#' 
#' @references
#' Shlomo, N., Antal, L., & Elliot, M. (2015). 
#' Measuring Disclosure Risk and Data Utility for Flexible Table Generators, 
#' Journal of Official Statistics, 31(2), 305-324. doi: \url{https://doi.org/10.1515/jos-2015-0019}
#'
#' @param f Vector of original counts
#' @param g Vector of perturbed counts
#'
#' @return  Hellinger distance or related utility measure
#' @export
#' 
#'
#' @examples
#' f <- 1:6
#' g <- c(0, 3, 3, 3, 6, 6)
#' print(c(
#'   HD = HD(f, g), 
#'   HDutility = HDutility(f, g), 
#'   maxdiff = max(abs(g - f)), 
#'   meanAbsDiff = mean(abs(g - f)), 
#'   rootMeanSquare = sqrt(mean((g - f)^2))
#' ))
HDutility <- function(f, g){ 
  1 - HD(f, g)/sqrt(sum(f))
}





