#' PLSrounding on portions of data at a time
#'
#' @param data data
#' @param loopId loopId 
#' @param ... dots
#' @param zeroCandidates zeroCandidates 
#' @param preRounded preRounded 
#' @param plsWeights plsWeights 
#' @param printInc printInc 
#' @param preDifference preDifference 
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
                            preRounded = NULL, 
                            plsWeights = NULL, 
                            printInc = TRUE,
                            preDifference = TRUE) {
  
  id <- unique(data[[loopId]])
  
  if (is.logical(preDifference)) {
    updatePreDifference <- preDifference
    preDifference <- NULL
  } else {
    stop("Supplied preDifference not implemented")  # Without stop Supplied preDifference will only be used when i=1
    updatePreDifference <- TRUE
  }
  
  a <- NULL
  a$publish <- preDifference
  
  for (i in seq_along(id)) {
    cat(sprintf("%4d: ", i))
    ai <- PLSrounding(data = data[data[[loopId]] == id[i], , drop = FALSE], 
                      ..., 
                      zeroCandidates = zeroCandidates, 
                      preRounded = preRounded, 
                      plsWeights = plsWeights,
                      printInc = printInc,
                      preDifference = if (updatePreDifference) a$publish else NULL)  # ifelse(updatePreDifference, a$publish, NULL) does not work since NULL
    if (i == 1) {
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





