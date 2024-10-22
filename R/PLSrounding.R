#' PLS inspired rounding
#' 
#'  Small count rounding of necessary inner cells are performed so that all small frequencies of cross-classifications to be published 
#' (publishable cells) are rounded. The publishable cells can be defined from a model formula, hierarchies or automatically from data.
#' 
#' This function is a user-friendly wrapper for \code{RoundViaDummy} with data frame output and with computed summary of the results.
#' See \code{\link{RoundViaDummy}} for more details.
#'
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (inner cells frequencies).  When \code{NULL} (default), microdata is assumed.
#' @param roundBase Rounding base
#' @param hierarchies List of hierarchies
#' @param formula Model formula defining publishable cells
#' @param dimVar The main dimensional variables and additional aggregating variables. This parameter can be  useful when hierarchies and formula are unspecified. 
#' @param maxRound Inner cells contributing to original publishable cells equal to or less than maxRound will be rounded
#' @param printInc Printing iteration information to console when TRUE  
#' @param output Possible non-NULL values are \code{"input"}, \code{"inner"} and \code{"publish"}. Then a single data frame is returned.
#' @param extend0  When `extend0` is set to `TRUE`, the data is automatically extended. 
#'                 This is relevant when `zeroCandidates = TRUE` (see \code{\link{RoundViaDummy}}).  
#'        Additionally, `extend0` can be specified as a list, representing the `varGroups` parameter 
#'        in the \code{\link[SSBtools]{Extend0}} function. 
#'        Can also be set to `"all"` which means that input codes in hierarchies are considered in addition to those in data.   
#' @param preAggregate When \code{TRUE}, the data will be aggregated beforehand within the function by the dimensional variables. 
#' @param aggregatePackage Package used to preAggregate. 
#'                         Parameter `pkg` to \code{\link[SSBtools]{aggregate_by_pkg}}.
#' @param aggregateNA Whether to include NAs in the grouping variables while preAggregate. 
#'                    Parameter `include_na` to \code{\link[SSBtools]{aggregate_by_pkg}}.
#' @param aggregateBaseOrder Parameter `base_order` to \code{\link[SSBtools]{aggregate_by_pkg}}, used when preAggregate.  
#'                           The default is set to `FALSE` to avoid unnecessary sorting operations.  
#'                           When `TRUE`, an attempt is made to return the same result with `data.table` as with base R.
#'                           This cannot be guaranteed due to potential variations in sorting behavior across different systems.
#' @param rowGroupsPackage Parameter `pkg` to \code{\link[SSBtools]{RowGroups}}.
#'               The parameter is input to \code{\link[SSBtools]{Formula2ModelMatrix}} 
#'               via \code{\link[SSBtools]{ModelMatrix}}. 
#' @param ... Further parameters sent to \code{RoundViaDummy}  
#'
#' @return Output is a four-element list with class attribute "PLSrounded", 
#'         which ensures informative printing and enables the use of \code{\link[SSBtools]{FormulaSelection}} on this object.
#'    \item{inner}{Data frame corresponding to input data with the main dimensional variables and with cell 
#'                frequencies (original, rounded, difference).}
#'    \item{publish}{Data frame of publishable data with the main dimensional variables and with cell frequencies 
#'                   (original, rounded, difference).}
#'    \item{metrics}{A named character vector of various statistics calculated from the two output data frames 
#'    ("\code{inner_}" used to distinguish). See examples below and the function \code{\link{HDutility}}.}
#'    \item{freqTable}{Matrix of frequencies of cell frequencies and absolute differences.
#'    For example, row "\code{rounded}" and column "\code{inn.4+}" is the number of rounded 
#'    inner cell frequencies greater than or equal to \code{4}.}
#'    
#' @seealso   \code{\link{RoundViaDummy}}, \code{\link{PLS2way}}, \code{\link[SSBtools]{ModelMatrix}} 
#' 
#' @references 
#' Langsrud, Ø. and Heldal, J. (2018): \dQuote{An Algorithm for Small Count Rounding of Tabular Data}. 
#' Presented at: \emph{Privacy in statistical databases}, Valencia, Spain. September 26-28, 2018.
#' \url{https://www.researchgate.net/publication/327768398_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data}
#' 
#' @encoding UTF8
#' 
#' @importFrom SSBtools CharacterDataFrame aggregate_by_pkg NamesFromModelMatrixInput 
#' @importFrom stats as.formula delete.response terms
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
#' # Using FormulaSelection()
#' FormulaSelection(a$publish, ~eu + year)
#' FormulaSelection(a, ~eu + year) # same as above
#' FormulaSelection(a)             # just a$publish
#' 
#' # Recalculation of maxdiff, HDutility, meanAbsDiff and rootMeanSquare
#' max(abs(a$publish[, "difference"]))
#' HDutility(a$publish[, "original"], a$publish[, "rounded"])
#' mean(abs(a$publish[, "difference"]))
#' sqrt(mean((a$publish[, "difference"])^2))
#' 
#' # Six lines below produce equivalent results 
#' # Ordering of rows can be different
#' PLSrounding(z, "freq") # All variables except "freq" as dimVar  
#' PLSrounding(z, "freq", dimVar = c("geo", "eu", "year"))
#' PLSrounding(z, "freq", formula = ~eu * year + geo * year)
#' PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eHrc"))
#' PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"))
#' PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo * year)
#' 
#' # Define publishable cells differently by making use of formula interface
#' PLSrounding(z, "freq", formula = ~eu * year + geo)
#' 
#' # Define publishable cells differently by making use of hierarchy interface
#' eHrc2 <- list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019"))
#' PLSrounding(z, "freq", hierarchies = eHrc2)
#' 
#' # Also possible to combine hierarchies and formula
#' PLSrounding(z, "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo + year)
#' 
#' # Single data frame output
#' PLSroundingInner(z, "freq", roundBase = 5, formula = ~geo + eu + year)
#' PLSroundingPublish(z, roundBase = 5, formula = ~geo + eu + year)
#' 
#' # Microdata input
#' PLSroundingInner(rbind(z, z), roundBase = 5, formula = ~geo + eu + year)
#' 
#' # Zero perturbed due to both  extend0 = TRUE and zeroCandidates = TRUE 
#' set.seed(12345)
#' PLSroundingInner(z[sample.int(5, 12, replace = TRUE), 1:3], 
#'                  formula = ~geo + eu + year, roundBase = 5, 
#'                  extend0 = TRUE, zeroCandidates = TRUE, printInc = TRUE)
#' 
#' # Parameter avoidHierarchical (see RoundViaDummy and ModelMatrix) 
#' PLSroundingPublish(z, roundBase = 5, formula = ~geo + eu + year, avoidHierarchical = TRUE)
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
#' 
#' # Use PLS2way to produce tables as in Langsrud and Heldal (2018) and to demonstrate 
#' # parameters maxRound, zeroCandidates and identifyNew (see RoundViaDummy).   
#' # Parameter rndSeed used to ensure same output as in reference.
#' exPSD <- SmallCountData("exPSD")
#' a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, rndSeed=124)
#' PLS2way(a, "original")  # Table 1
#' PLS2way(a)  # Table 2
#' a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, identifyNew = FALSE, rndSeed=124)
#' PLS2way(a)  # Table 3
#' a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, maxRound = 7)
#' PLS2way(a)  # Values in col1 rounded
#' a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, zeroCandidates = TRUE)
#' PLS2way(a)  # (row3, col4): original is 0 and rounded is 5
PLSrounding <- function(data, freqVar = NULL, roundBase = 3, hierarchies = NULL, formula = NULL, 
                        dimVar = NULL,
                        maxRound = roundBase-1, printInc = nrow(data)>1000, 
                        output = NULL, 
                        extend0 = FALSE,
                        preAggregate = is.null(freqVar),
                        aggregatePackage = "base",
                        aggregateNA = TRUE,
                        aggregateBaseOrder = FALSE,
                        rowGroupsPackage = aggregatePackage,
                        ...) {
  
  
  force(preAggregate)
  
  names_data <- names(data)
  
  if(!is.null(output)){
    if(!(output %in% c("input", "inner", "publish")))
      stop('Allowed non-NULL values of parameter output are "input", "inner" and "publish".')
    
  } else {
    output <- ""
  }
  
  isExtend0 <- IsExtend0(extend0)
  
  if (preAggregate | output == "input" | isExtend0) {
    if (printInc & preAggregate) {
      cat("[preAggregate ", dim(data)[1], "*", dim(data)[2], "->", sep = "")
      flush.console()
    }
    dVar <- NamesFromModelMatrixInput(hierarchies = hierarchies, formula = formula, dimVar = dimVar)
    if (!length(dVar)) {
          if (is.null(freqVar)){
            dVar <- names(data)
          } else {
            freqVarName <- names(data[1, freqVar, drop = FALSE])
            dVar <- names(data[1, !(names(data) %in% freqVarName), drop = FALSE])
          }
    }
    dVar <- unique(dVar)
    
    if (preAggregate) {
      if (is.null(freqVar)) {
        #data <- aggregate(list(f_Re_qVa_r = data[[dVar[1]]]), data[dVar], length)
        data <- aggregate_by_pkg(
          data = data,
          by = dVar,
          var = dVar[1],
          pkg =  aggregatePackage,
          include_na = aggregateNA,
          fun = length,
          base_order = aggregateBaseOrder)
        freqVar <- "f_Re_qVa_r"
        names(data)[length(dVar) + 1] <- freqVar 
      } else {
        #data <- aggregate(data[freqVar], data[dVar], sum)
        data <- aggregate_by_pkg(
          data = data,
          by = dVar,
          var = freqVar,
          pkg =  aggregatePackage,
          include_na = aggregateNA,
          fun = sum,
          base_order = aggregateBaseOrder)
      }
      if (printInc) {
        cat(dim(data)[1], "*", dim(data)[2], "]\n", sep = "")
        flush.console()
      }
    }
  }
  
  if (isExtend0) {
    if (printInc) {
      cat("[extend0 ", dim(data)[1], "*", dim(data)[2], "->", sep = "")
      flush.console()
    }
    
    data <- Extend0fromModelMatrixInput(data = data, 
                                        freqName = freqVar, 
                                        hierarchies = hierarchies, 
                                        formula = formula,
                                        dimVar = dimVar,
                                        extend0 = extend0, 
                                        dVar = dVar, ...)
    
    
    if (printInc) {
      cat(dim(data)[1], "*", dim(data)[2], "]\n", sep = "")
      flush.console()
    }
  }
  
  if (output == "input"){
    if (!is.null(freqVar)) {
      if (freqVar == "f_Re_qVa_r") {
        if (!("freq" %in% names_data)) {
          names(data)[names(data) == freqVar] <- "freq"
          freqVar <- "freq"
        }
      }
    }
    inner <- list(inner = data)
    return(data[, c(dVar, freqVar), drop = FALSE])
  }
  
  
  if(!is.null(list(...)$Version)){   # For testing
    z <- RoundViaDummy_Version_0.3.0(data = data, freqVar = freqVar, formula = formula, roundBase = roundBase, hierarchies = hierarchies, ...) 
  } else {
    z <- RoundViaDummy(data = data, freqVar = freqVar, formula = formula, roundBase = roundBase, hierarchies = hierarchies,
                       dimVar = dimVar,
                     maxRound = maxRound, printInc = printInc,
                     rowGroupsPackage = rowGroupsPackage, ...)
  }
  
  
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
  
  freqTable <- cbind(TabCat(z$yInner, roundBase, "inn.", maxRound), TabCat(z$yPublish, roundBase, "pub.", maxRound))
  
  # freqVarName <- names(data[1, freqVar, drop = FALSE])
  
  out <- NULL
  
  if (!is.null(z$crossTable)) {
    cNames <- colnames(data)[colnames(data) %in% colnames(z$crossTable)]
    if (output != "publish"){
      out$inner <- cbind( CharacterDataFrame(data[, cNames, drop = FALSE]), z$yInner)
    }
    if (output != "inner"){
      out$publish <- cbind(as.data.frame(z$crossTable[, cNames, drop = FALSE], stringsAsFactors = FALSE), z$yPublish)
      rownames(out$publish) <- NULL
      startRow <- attr(z$crossTable, "startRow", exact = TRUE)
      if (!is.null(startRow)) {
        attr(out$publish, "startRow") <- startRow
      }
    }
  } else {
    if (output != "publish"){
      out$inner <- as.data.frame(z$yInner)
    }
    if (output != "inner"){
      out$publish <- as.data.frame(z$yPublish)
    }
  }
  
  if (output == "inner"){
    return(out$inner)
  }
  
  if (output == "publish"){
    return(out$publish)
  }
  
  
  if (any(duplicated(colnames(out$inner)))) {
    warning(paste("Duplicated colnames in output:", paste(colnames(out$inner)[duplicated(colnames(out$inner))], collapse = ", ")))
  }
  
  
  out$metrics <- c(roundBase = roundBase, maxRound = maxRound, maxdiff = maxdiff, 
                   inner_HDutility = inner_HDutility, HDutility = publish_HDutility, 
                   inner_meanAbsDiff = inner_meanAbsDiff, meanAbsDiff = publish_meanAbsDiff, 
                   inner_rootMeanSquare = inner_rootMeanSquare, rootMeanSquare = publish_rootMeanSquare)
  out$freqTable <- freqTable
  
  if (!is.null(z$x))
    out$x <- z$x
    
  out
  return(structure(out, class = "PLSrounded"))
}



#' @rdname PLSrounding
#' @export
PLSroundingInner  <- function(..., output = "inner") {
  PLSrounding(..., output = output)
}


#' @rdname PLSrounding
#' @export
PLSroundingPublish  <- function(..., output = "publish") {
  PLSrounding(..., output = output)
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


#' FormulaSelection  method for PLSrounded
#'
#' @param x PLSrounded object 
#' @param formula `formula` parameter to \code{\link[SSBtools]{FormulaSelection}}.
#'        When `NULL` (default), the publish data frame is returned without any limitation.  
#' @param intercept `intercept` parameter to `FormulaSelection`.
#' @param logical `logical` parameter to `FormulaSelection`.
#'
#' @return Limited version of the publish data frame
#' @importFrom SSBtools FormulaSelection
#' @export
FormulaSelection.PLSrounded <- function(x, formula = NULL, intercept = NA, logical = FALSE) {
  if (is.null(formula)) {
    return(x$publish)
  }
  SSBtools::FormulaSelection(x$publish, formula, intercept, logical)
}


TabCat <- function(z, b, s, m) {
  x <- rbind(original = table(ToCat(z[, "original"], b, m)), 
             rounded = table(ToCat(z[, "rounded"], b, m)), 
             absDiff = table(ToCat(abs(z[, "difference"]), b, m)))
  x <- x[, colnames(x) != "0.5", drop = FALSE]
  x <- cbind(x, all = NROW(z))
  colnames(x) <- paste(s, colnames(x), sep = "")
  x
}


ToCat <- function(x, b, m) {
  if(m >= b)
    return(ToCatm(x, b, m))
  x[x > b] <- b + 1
  x[x > 0 & x < b] <- 1 - as.numeric(b == 1)/2
  x <- factor(x, levels = c(0, 1 - as.numeric(b == 1)/2, b, b + 1))
  levels(x)[4] <- paste(levels(x)[4], "+", sep = "")
  if (b > 2) 
    levels(x)[2] <- paste(1, as.character(b - 1), sep = "-")
  x
}


ToCatm <- function(x, b, m) {
  x[x > m] <- m + 1
  x[x > (b-1) & x<=m ] <- b 
  x[x > 0 & x < b] <- 1 - as.numeric(b == 1)/2
  x <- factor(x, levels = c(0, 1 - as.numeric(b == 1)/2, b, m + 1))
  levels(x)[4] <- paste(levels(x)[4], "+", sep = "")
  if (b > 2) 
    levels(x)[2] <- paste(1, as.character(b - 1), sep = "-")
  if(m > b)
    levels(x)[3] <- paste(as.character(b), as.character(m), sep = "-")
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
#' Journal of Official Statistics, 31(2), 305-324. \doi{10.1515/jos-2015-0019}
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





