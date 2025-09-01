#' PLSrounding on portions of data at a time
#' 
#' The \code{\link{PLSrounding}} runs are coordinated by using preliminary differences as input for the next run (parameter `preDifference`)
#' 
#' Note that in this function `zeroCandidates`, `forceInner`, `preRounded` and `plsWeights` cannot be supplied as vectors.
#' They may be specified as functions or as variables in the input data.
#'
#' @param data Input data as a data frame (inner cells)
#' @param loopId Variable holding id for loops
#' @param ... `PLSrounding` parameters 
#' @param zeroCandidates `PLSrounding` parameter (see details) 
#' @param forceInner `PLSrounding` parameter (see details)
#' @param preRounded `PLSrounding` parameter (see details)
#' @param plsWeights `PLSrounding` parameter (see details)
#' @param printInc Printing iteration information to console when TRUE
#' @param preDifference When TRUE, the `preDifference` parameter to `PLSrounding` is used. Each time with the differences obtained so far.  
#' @param preOutput preOutput The function can continue from output from a previous run
#' @param rndSeed If non-NULL, a random generator seed to be set locally at the beginning of `PLSroundingLoop` without affecting the random value stream in R.
#'                Within `PLSroundingLoop`, `PLSrounding` is called with `rndSeed = NULL`.
#' @param action_unused_dots `PLSrounding` parameter.             
#'
#' @return As output from \code{\link{PLSrounding}} 
#' @export
#'
#' @examples
#' mf2 <- ~region + fylke * hovedint
#' z2 <- SmallCountData("z2")
#' a <- PLSroundingLoop(z2, loopId = "kostragr", freqVar = "ant", formula = mf2)
#' a
PLSroundingLoop <- function(data, 
                            loopId, 
                            ..., 
                            zeroCandidates = FALSE, 
                            forceInner = FALSE,
                            preRounded = NULL, 
                            plsWeights = NULL, 
                            printInc = TRUE,
                            preDifference = TRUE,
                            preOutput = NULL,
                            rndSeed = 123,
                            action_unused_dots = "warn"
                            ) {
  if (!is.null(rndSeed)) {
    if (!exists(".Random.seed")) 
      if (runif(1) < 0) 
        stop("Now seed exists")
    exitSeed <- .Random.seed
    on.exit(.Random.seed <<- exitSeed)
    set.seed(rndSeed)
  }
  
  id <- unique(data[[loopId]])
  
  if (is.logical(preDifference)) {
    updatePreDifference <- preDifference
    preDifference <- NULL
  } else {
    stop("Supplied preDifference not implemented. Use preOutput instead.")  # Without stop Supplied preDifference will only be used when i=1
    updatePreDifference <- TRUE
  }
  
  preOutputInInput <- as.numeric(!is.null(preOutput))
  
  a <- preOutput
  
  lengthStop <- FALSE
  
  if (length(zeroCandidates) > 1) lengthStop <- TRUE
  if (length(forceInner) > 1) lengthStop <- TRUE
  if (length(preRounded) > 1) lengthStop <- TRUE
  if (length(plsWeights) > 1) lengthStop <- TRUE
  
  if (lengthStop) stop("zeroCandidates, forceInner, preRounded and plsWeights cannot be supplied as vectors in PLSroundingLoop.")
  
  
  for (i in seq_along(id)) {
    if (i == 2) {
      action_unused_dots <- "none"
    }
    if(printInc) cat(sprintf("%4d: ", i))
    ai <- PLSrounding(data = data[data[[loopId]] == id[i], , drop = FALSE], 
                      ..., 
                      zeroCandidates = zeroCandidates, 
                      forceInner = forceInner,
                      preRounded = preRounded, 
                      plsWeights = plsWeights,
                      printInc = printInc,
                      preDifference = if (updatePreDifference) a$publish else NULL,  # ifelse(updatePreDifference, a$publish, NULL) does not work since NULL
                      rndSeed = NULL, 
                      action_unused_dots = action_unused_dots) 
    if (i == (1 - preOutputInInput)) {
      a <- ai
    } else {
      a <- UpdatePLSrounded(a, ai)
    }
  }
  a
}


UpdatePLSrounded <- function(a, b) {
  if (!identical(colnames(a$publish), colnames(b$publish))) {
    stop("colnames mismatch")
  }
  a$inner <- CombinePLSrounded(a$inner, b$inner)
  a$publish <- CombinePLSrounded(a$publish, b$publish)
  metricsfreqTable <- UpdateMetricsfreqTable(a)
  a$metrics <- metricsfreqTable$metrics
  a$freqTable <- metricsfreqTable$freqTable
  if (!is.null(a$x)) {
    warning("x not updated when looping (xReturn)")
  }
  a
}

CombinePLSrounded <- function(a, b) {
  dimVar <- colnames(a)
  numVar <- c("original", "rounded", "difference")
  dimVar <- dimVar[!(dimVar %in% numVar)]
  
  ma <- Match(b[dimVar], a[dimVar])
  
  a[ma[!is.na(ma)], numVar] <- a[ma[!is.na(ma)], numVar, drop = FALSE] + b[!is.na(ma), numVar, drop = FALSE]
  
  rbind(a, b[is.na(ma), , drop = FALSE])
}

UpdateMetricsfreqTable <- function(a) {
  roundBase <- as.numeric(a$metrics["roundBase"])  # naming 
  maxRound <- as.numeric(a$metrics["maxRound"])
  maxdiff <- max(abs(a$publish[, "difference"]))
  
  inner_HDutility <- HDutility(a$inner[, "original"], a$inner[, "rounded"])
  publish_HDutility <- HDutility(a$publish[, "original"], a$publish[, "rounded"])
  inner_meanAbsDiff <- mean(abs(a$inner[, "difference"]))
  publish_meanAbsDiff <- mean(abs(a$publish[, "difference"]))
  inner_rootMeanSquare <- sqrt(mean((a$inner[, "difference"])^2))
  publish_rootMeanSquare <- sqrt(mean((a$publish[, "difference"])^2))
  freqTable <- cbind(TabCat(a$inner, roundBase, "inn.", maxRound), TabCat(a$publish, roundBase, "pub.", maxRound))
  out <- NULL
  out$metrics <- c(roundBase = roundBase, maxRound = maxRound, maxdiff = maxdiff, 
                   inner_HDutility = inner_HDutility, HDutility = publish_HDutility, 
                   inner_meanAbsDiff = inner_meanAbsDiff, meanAbsDiff = publish_meanAbsDiff, 
                   inner_rootMeanSquare = inner_rootMeanSquare, 
                   rootMeanSquare = publish_rootMeanSquare)
  out$freqTable <- freqTable
  out
}





