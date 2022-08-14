#' Small Count Rounding of Tabular Data
#'
#' Small count rounding via a dummy matrix and by an algorithm inspired by PLS
#'
#' Small count rounding of necessary inner cells are performed so that all small frequencies of cross-classifications to be published 
#' (publishable cells) are rounded. This is equivalent to changing micro data since frequencies of unique combinations are changed. 
#' Thus, additivity and consistency are guaranteed. The matrix multiplication formula is: 
#' \code{yPublish} \code{=} \code{t(x)} \code{\%*\%}  \code{yInner}, where \code{x} is the dummy matrix. 
#' 
#' 
#' Parameters `zeroCandidates`, `forceInner`,  `preRounded` and `plsWeights` can be specified as functions.
#' The supplied functions take the following arguments: `data`, `yPublish`,  `yInner`, `crossTable`, `x`, `roundBase`, `maxRound`, and `...`, 
#'                   where the first two are numeric vectors of original counts. 
#' When `allSmall` is `TRUE`,  `forceInner` is set to  `function(yInner, maxRound, ...)` `yInner <= maxRound`.                   
#'
#' @encoding UTF8
#' @md
#'
#' @param data Input data as a data frame (inner cells)
#' @param freqVar Variable holding counts (name or number)
#' @param formula Model formula defining publishable cells. Will be used to calculate \code{x} (via \code{\link{ModelMatrix}}). 
#' When NULL, x must be supplied.
#' @param roundBase Rounding base
#' @param singleRandom Single random draw when TRUE (instead of algorithm)
#' @param crossTable	When TRUE, cross table in output and caculations via FormulaSums()
#' @param total	String used to name totals
#' @param maxIterRows See details
#' @param maxIter Maximum number of iterations
#' @param x Dummy matrix defining publishable cells
#' @param hierarchies List of hierarchies, which can be converted by \code{\link{AutoHierarchies}}. 
#'        Thus, a single string as hierarchy input is assumed to be a total code. 
#'        Exceptions are \code{"rowFactor"} or \code{""}, which correspond to only using the categories in the data.
#' @param xReturn Dummy matrix in output when TRUE (as input parameter \code{x})
#' @param maxRound Inner cells contributing to original publishable cells equal to or less than maxRound will be rounded.
#' @param zeroCandidates When TRUE, inner cells in input with zero count (and multiple of roundBase when maxRound is in use)
#'             contributing to publishable cells will be included as candidates to obtain roundBase value. 
#'             With vector input, the rule is specified individually for each cell. 
#'             This can be specified as a vector, a variable in data or a function generating it (see details).   
#' @param forceInner When TRUE, all inner cells will be rounded. Use vector input to force individual cells to be rounded.
#'                   This can be specified as a vector, a variable in data or a function generating it (see details).
#'                   Can be combined with parameter zeroCandidates to allow zeros and roundBase multiples to be rounded up.
#' @param identifyNew  When TRUE, new cells may be identified after initial rounding to ensure that no nonzero rounded 
#'        publishable cells are less than roundBase.  
#' @param step When \code{step>1}, the original forward part of the algorithm is replaced by a kind of stepwise. 
#'       After \code{step} steps forward, backward steps may be performed. The \code{step} parameter is also used 
#'       for backward-forward iteration at the end of the algorithm; \code{step} backward steps may be performed.
#' @param preRounded A vector or a variable in data that contains a mixture of missing values and predetermined values of rounded inner cells. 
#'                   Can also be specified as a function generating it (see details).
#' @param leverageCheck When TRUE, all inner cells that depends linearly on the published cells and with small frequencies
#'        (\code{<=maxRound}) will be rounded. 
#'        The computation of leverages can be very time and memory consuming. 
#'        The function \code{\link{Reduce0exact}} is called. 
#'        The default leverage limit is `0.999999`. Another limit can be sent as input instead of `TRUE`.  
#'        Checking is performed before and after (since new zeros) rounding. Extra iterations are performed when needed.  
#' @param easyCheck A light version of the above leverage checking. 
#'                  Checking is performed after rounding. Extra iterations are performed when needed.
#'                  `Reduce0exact` is called with `reduceByLeverage=FALSE` and `reduceByColSums=TRUE`.
#' @param printInc Printing iteration information to console when TRUE
#' @param rndSeed If non-NULL, a random generator seed to be used locally within the function without affecting the random value stream in R.      
#' @param dimVar The main dimensional variables and additional aggregating variables. This parameter can be  useful when hierarchies and formula are unspecified.   
#' @param plsWeights A vector of weights for each cell to be published or a function generating it (see details). For use in the algorithm criterion.
#' @param preDifference   A data.frame with differences already obtained from rounding another subset of data. 
#'                        There must be columns that match `crossTable`. Differences must be in the last column.
#' @param allSmall When TRUE, all small inner cells (`<= maxRound`) are rounded. This parameter is a simplified alternative to specifying `forceInner`  (see details).                            
#' @param ... Further parameters sent to \code{\link{ModelMatrix}}.
#'            In particular, one can specify `removeEmpty=TRUE` to omit empty combinations.     
#'            The parameter `inputInOutput` can be used to specify whether to include codes from input.
#'            The parameter `avoidHierarchical` (\code{\link{Formula2ModelMatrix}}) can be combined with formula input.
#' @note Iterations are needed since after initial rounding of identified cells, new cells are identified.
#' If cases of a high number of identified cells the algorithm can be too memory consuming (unless singleRandom=TRUE).
#' To avoid problems, not more than maxIterRows cells are rounded in each iteration.
#' The iteration limit (maxIter) is by default set to be high since a low number of maxIterRows may need a high number of iterations.
#'
#' @return A list where the two first elements are two column matrices.
#' The first matrix consists of inner cells and the second of cells to be published.
#' In each matrix the first and the second column contains, respectively, original and rounded values.
#' By default the cross table is the third element of the output list.
#' 
#' @seealso  See the  user-friendly wrapper \code{\link{PLSrounding}}
#'   and see \code{Round2} for rounding by other algorithm
#' @importFrom stats as.formula hat model.frame model.matrix
#' @importFrom SSBtools FormulaSums matlabColon Hierarchies2ModelMatrix FindHierarchies Reduce0exact MakeFreq  ModelMatrix As_TsparseMatrix
#' @importFrom utils flush.console
#' @importFrom  Matrix Matrix sparse.model.matrix Diagonal
#' @importFrom  methods hasArg
#' @export
#'
#' @examples
#' # See similar and related examples in PLSrounding documentation
#' RoundViaDummy(SmallCountData("e6"), "freq")
#' RoundViaDummy(SmallCountData("e6"), "freq", formula = ~eu * year + geo)
#' RoundViaDummy(SmallCountData("e6"), "freq", hierarchies = 
#'    list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019")))
#' 
#' RoundViaDummy(SmallCountData('z2'), 
#'               'ant', ~region + hovedint + fylke*hovedint + kostragr*hovedint, 10)
#' mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
#' a <- RoundViaDummy(SmallCountData('z3'), 'ant', mf, 5)
#' b <- RoundViaDummy(SmallCountData('sosialFiktiv'), 'ant', mf, 4)
#' print(cor(b[[2]]),digits=12) # Correlation between original and rounded
#' 
#' # Demonstrate parameter leverageCheck 
#' # The 42nd inner cell must be rounded since it can be revealed from the published cells.
#' mf2 <- ~region + hovedint + fylke * hovedint + kostragr * hovedint
#' RoundViaDummy(SmallCountData("z2"), "ant", mf2, leverageCheck = FALSE)$yInner[42, ]
#' RoundViaDummy(SmallCountData("z2"), "ant", mf2, leverageCheck = TRUE)$yInner[42, ]
#' 
#' \dontrun{
#' # Demonstrate parameters maxRound, zeroCandidates and forceInner 
#' # by tabulating the inner cells that have been changed.
#' z4 <- SmallCountData("sosialFiktiv")
#' for (forceInner in c("FALSE", "z4$ant < 10")) 
#'   for (zeroCandidates in c(FALSE, TRUE)) 
#'     for (maxRound in c(2, 5)) {
#'       set.seed(123)
#'       a <- RoundViaDummy(z4, "ant", formula = mf, maxRound = maxRound, 
#'                          zeroCandidates = zeroCandidates, 
#'                          forceInner = eval(parse(text = forceInner)))
#'       change <- a$yInner[, "original"] != a$yInner[, "rounded"]
#'       cat("\n\n---------------------------------------------------\n")
#'       cat("      maxRound:", maxRound, "\n")
#'       cat("zeroCandidates:", zeroCandidates, "\n")
#'       cat("    forceInner:", forceInner, "\n\n")
#'       print(table(original = a$yInner[change, "original"], rounded = a$yInner[change, "rounded"]))
#'       cat("---------------------------------------------------\n")
#'     }
#' }
RoundViaDummy <- function(data, freqVar, formula = NULL, roundBase = 3, singleRandom = FALSE,
                          crossTable=TRUE, total = "Total",  maxIterRows = 1000, maxIter = 1E7,
                          x = NULL,  hierarchies = NULL, xReturn = FALSE, maxRound = roundBase-1,
                          zeroCandidates = FALSE, forceInner = FALSE, identifyNew = TRUE, step = 0,
                          preRounded = NULL,
                          leverageCheck = FALSE, 
                          easyCheck = TRUE,
                          printInc = TRUE,  rndSeed = 123, dimVar = NULL, 
                          plsWeights = NULL, 
                          preDifference = NULL, 
                          allSmall = FALSE, 
                          ...) {
  
  if (!is.null(rndSeed)) {
    if (!exists(".Random.seed")) 
      if (runif(1) < 0) 
        stop("Now seed exists")
    exitSeed <- .Random.seed
    on.exit(.Random.seed <<- exitSeed)
    set.seed(rndSeed)
  }
  
  if (hasArg("roundbase"))
    stop('Misspelled parameter "roundbase" found. Use "roundBase".')
  
  if(roundBase<1){
    stop(paste("roundBase =", roundBase,"is not allowed"))
  }
  
  if(roundBase == 1L){
    warning("Special algorithm when roundBase is 1, forceInner and zeroCandidates set to TRUE.")
    zeroCandidates = TRUE
    forceInner = TRUE 
    identifyNew = FALSE
  }
  
  maxBase <- maxRound + 1
  
  if(identifyNew)
    if(maxBase<roundBase)
      stop("maxRound cannot be smaller than roundBase-1 when identifyNew is TRUE")

  if (is.logical(leverageCheck)) {
    leverageCheck <- 0.999999 * as.numeric(leverageCheck)
  }
  
  if(printInc){
    cat("[")
    flush.console()
  }
  
  
  if (!is.null(formula) & is.null(hierarchies) & is.null(x)){
    
    if(length(total)>1){
      total <- total[1]
      warning("Only first element of total is used when formula input.")
    }
  }
  
  
  if(is.null(freqVar)){       # microdata (freq=1) when empty freqVar 
    freqVar <- "f_Re_qVa_r"
    data[[freqVar]] <- 1L
  }
  
  
  if (is.null(formula) & is.null(hierarchies) & is.null(x) & is.null(dimVar)){
    freqVarName <- names(data[1, freqVar, drop = FALSE])
    dimVar <- names(data[1, !(names(data) %in% freqVarName), drop = FALSE])
  }
  
  
  x = ModelMatrix(data = data, hierarchies=hierarchies, formula =formula, crossTable=crossTable, modelMatrix = x, total=total, dimVar=dimVar, ...)
  
  #if(is.null(freqVar) & xReturn){
  #  return(x)
  #}
  
  if(is.list(x)){
    crossTable = TRUE
    crossTab <- x$crossTable
    x <- x$modelMatrix
  } else {
    crossTable = FALSE
    crossTab <- NULL
  }
  
  
  yInner <- data[, freqVar, drop = TRUE]
  
  
  yPublish <- Matrix::crossprod(x, yInner)[, 1, drop = TRUE]
  
  

  ########### zeroCandidates ###########
  if (is.function(zeroCandidates)) {
    zeroCandidates <- zeroCandidates(data = data, yPublish = yPublish, yInner = yInner, crossTable = crossTab, x = x, roundBase = roundBase, maxRound = maxRound, ...)
  }
  if (length(zeroCandidates) == 1 & !is.logical(zeroCandidates[1])) {
    zeroCandidates <- data[[zeroCandidates]]
  } else {
    if (length(zeroCandidates) != 1 & length(zeroCandidates) != length(yInner)) {
      stop("Wrong length of zeroCandidates")
    }
  }
  

  ########### forceInner ###########
  if (allSmall) {
    if (!identical(forceInner, FALSE)) {
      warning("forceInner ignored since allSmall is TRUE")
    }
    forceInner <- function(yInner, maxRound, ...) yInner <= maxRound
  }
  
  if (is.function(forceInner)) {
    forceInner <- forceInner(data = data, yPublish = yPublish, yInner = yInner, crossTable = crossTab, x = x, roundBase = roundBase, maxRound = maxRound, ...)
  }
  if (length(forceInner) == 1 & !is.logical(forceInner[1])) {
    forceInner <- data[[forceInner]]
  } else {
    if (length(forceInner) != 1 & length(forceInner) != length(yInner)) {
      stop("Wrong length of forceInner")
    }
  }
  
  
  ########### preRounded ###########
  if (!is.null(preRounded)) {
    if (is.function(preRounded)) {
      preRounded <- preRounded(data = data, yPublish = yPublish, yInner = yInner, crossTable = crossTab, x = x, roundBase = roundBase, maxRound = maxRound, ...)
    }
    if (length(preRounded) == 1) {
      preRounded <- data[[preRounded]]
    } else {
      if (length(preRounded) != length(yInner)) {
        stop("Wrong length of preRounded")
      }
    }
  }
  
  
  ########### plsWeights ###########
  if(!is.null(plsWeights)){
    if (is.function(plsWeights)){
      plsWeights <- plsWeights(data = data, yPublish = yPublish, yInner = yInner, crossTable = crossTab, x = x, roundBase = roundBase, maxRound = maxRound, ...)
    }
    if(length(plsWeights)!= length(yPublish)){
      stop("Wrong length of plsWeights")
    }
  }
  
  
  
  
  if (!is.null(preDifference)) {
    preDifferences <- rep(0L, length(yPublish))
    ma <- Match(crossTab, preDifference[names(crossTab)])
    preDifferences[!is.na(ma)] <- preDifference[[ncol(preDifference)]][ma[!is.na(ma)]]
    
    maxCsum <- Matrix::colSums(x[, !is.na(ma), drop = FALSE])
    maxCsumInd <- which.max(maxCsum)
    maxCsum <- maxCsum[maxCsumInd]
    if (maxCsum != nrow(x)) {
      stop("Overall total needed in preDifference")
    }
    preDifferenceTot <- preDifferences[!is.na(ma)][maxCsumInd]
    if (is.null(preRounded) & abs(preDifferenceTot) > roundBase) {
      warning("Big overall total in preDifference. preRounded-input may trig a better algorithm")
    }
  } else {
    preDifferences <- NULL
    preDifferenceTot <- NULL
  }
  
  a <- PlsRoundSparse(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish, singleRandom = singleRandom,maxIter=maxIter, maxIterRows=maxIterRows, 
                      maxBase = maxBase, zeroCandidates = zeroCandidates, forceInner = forceInner, identifyNew = identifyNew, step = step, 
                      preRounded = preRounded,
                      easyCheck = easyCheck,
                      leverageCheck = leverageCheck,  printInc = printInc, plsWeights = plsWeights,
                      preDifferences = preDifferences, preDifferenceTot = preDifferenceTot)
  if(printInc){
    cat("]\n")
    flush.console()
  }
  
  if (xReturn){ # copy of code below with x as extra
    if(crossTable){
      return(list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
                  yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]),
                  crossTable = crossTab, x = x))
    } else {
      return(list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
         yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]), x = x))
    }
  }
  
  
  if(crossTable)
    return(list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
                yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]),
                crossTable = crossTab))
  
  list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
       yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]))
}



IntegerCbind = function(original,rounded){ # To ensure integer when integer input
  if(is.integer(original))
    return(cbind(original=original,rounded = as.integer(rounded)))
  cbind(original=original,rounded=rounded)
}




PlsRoundSparse <- function(x, roundBase = 3, yInner, yPublish = Matrix::crossprod(x, yInner)[, 1, drop = TRUE],
                           singleRandom = FALSE, maxIter = 1E6, maxIterRows = 1000, 
                           maxBase = roundBase, zeroCandidates = FALSE, forceInner = FALSE,
                           identifyNew = TRUE, step = 0,
                           forceFromFirstIter = TRUE, 
                           preRounded = NULL,
                           plsWeights = NULL,
                           easyCheck = TRUE, leverageCheck = 0, printInc = TRUE,
                           preDifferences = NULL, preDifferenceTot = NULL) { # maxIter henger sammen med maxIterRows

  

  if (!is.null(preDifferences)) {
    yInnerExact <- sum(yInner) - preDifferenceTot  # Possible since sum is only usage of yInnerExact
    yPublishExact <- yPublish - preDifferences
  } else {
    yInnerExact <- yInner
    yPublishExact <- yPublish
  }

  if(any(zeroCandidates)){
    if(length(zeroCandidates)==1){
      zeroCandidates =  (yInner %% roundBase) == 0
    } else {
      zeroCandidates <- zeroCandidates & (yInner %% roundBase) == 0
    }
  } else {
    zeroCandidates <- rep(FALSE, length(yInner))
  }    
  
  leverage <- NULL
  
  if (leverageCheck) {
    if (any(!forceInner) & any(yInner < maxBase) & leverageCheck!=0.99999901234) {
      if(printInc){
        cat("{")
        flush.console()
      }
      leverage <- Reduce0exact(x, z = Matrix(yPublish,ncol=1),y = Matrix(yInner,ncol=1), reduceByLeverage = TRUE, leverageLimit = leverageCheck, printInc = printInc)$yKnown
      
      leverage[!(yInner < maxBase & ( (yInner %% roundBase) > 0))] <- FALSE
      
      if(any(forceInner))
        leverage[forceInner] <- FALSE
      
      if(printInc){
        cat("}")
        flush.console()
      }
    } 
  }
  
  if (is.null(leverage))
    leverage <- rep(FALSE, length(yInner))
  
  if(any(forceInner)){
    if(length(forceInner)==1){
      forceInner <- rep(TRUE, length(yInner))
    } else {
      if(length(yInner) != length(forceInner))
        stop("Wrong length of forceInner")
    }
    forceInner[((yInner %% roundBase) == 0) & !zeroCandidates] <- FALSE
  } else {
    forceInner <- rep(FALSE, length(yInner))
  }
  

  if(forceFromFirstIter){
    create_supRowsForce <- TRUE
  } else {
    create_supRowsForce <- (maxBase != roundBase) |  !identifyNew  | any(zeroCandidates) | any(forceInner) | any(leverage)
  }
  
  
  if(create_supRowsForce){  
    maxBase = as.integer(maxBase)
    suppPublish = yPublish < maxBase & yPublish > 0
    suppInput <- yInner < maxBase & ( (yInner %% roundBase) > 0 | zeroCandidates )#  Indre celler med verdier som er 'undertrykkbare'
    supRowsForce <- (Matrix::rowSums(x[, suppPublish, drop = FALSE]) > 0 & suppInput) | forceInner
    
    if(any(!(supRowsForce[leverage]))){
      if(printInc) {
        cat(paste0("(check:",sum(!(supRowsForce[leverage])),")"))
      }
      supRowsForce[leverage] <- TRUE
    }
    
  } else {
    supRowsForce <- rep(FALSE, length(yInner))
  }
  
  
  if (!is.null(preRounded)) {
    yInner[!is.na(preRounded)] <- preRounded[!is.na(preRounded)]
    yPublish <- Matrix::crossprod(x, yInner)[, 1, drop = TRUE]
    preRounded <- !is.na(preRounded)
  } else {
    preRounded <- FALSE
  }

  if(any(preRounded))
    supRowsForce[preRounded] = FALSE
  

  i = 0
  while (i<maxIter) {
    i = i+1
    if (i == 1)
      a <- PlsRoundSparseSingle(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish,
                                singleRandom = singleRandom, yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows, 
                                supRowsForce = supRowsForce, identifyNew = !forceFromFirstIter, 
                                preRounded = preRounded,
                                step = step, printInc = printInc, plsWeights = plsWeights)
    else
      a <- PlsRoundSparseSingle(x = x, roundBase = roundBase, yInner = a[[1]], yPublish = yPublish,
                                singleRandom = singleRandom, 
                                suppPublish = suppRoundPublish,
                                yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows, 
                                supRowsForce = supRowsForce, identifyNew = identifyNew,
                                preRounded = preRounded,
                                step = step, printInc = printInc, plsWeights = plsWeights)
    
      yPublish  <- a[[2]][, 1, drop = TRUE]
      suppRoundPublish <- yPublish < roundBase & yPublish > 0
      supRowsForce[a[[3]]] <- FALSE
      
      return_a <- FALSE
      
      if(!identifyNew){
        if (!any(supRowsForce) )
          return_a <- TRUE
      } else {
        if (!any(suppRoundPublish) & !any(supRowsForce) ){
          return_a <- TRUE
        } else {
          if (sum(a[[3]])==0 & !any(supRowsForce) ){
            return_a <- TRUE
            if(any(preRounded))
              warning("Due to preRounded, not all published cells may be sufficiently rounded")
            else
              warning("Something wrong, ... not all published cells may be sufficiently rounded")
          }
        }
      }
      
      if (return_a) {
        if ( easyCheck | leverageCheck) { 
          if((printInc  & as.logical(leverageCheck))) cat("{")
          leverage <- Reduce0exact(x, z = Matrix(yPublish,ncol=1),y = Matrix(a[[1]],ncol=1), reduceByLeverage = as.logical(leverageCheck), 
                                   leverageLimit = leverageCheck, reduceByColSums=TRUE, 
                                   printInc =  (printInc  & as.logical(leverageCheck)))$yKnown
          if((printInc  & as.logical(leverageCheck))) {
            cat("}")
            flush.console()
          }
          if(any(leverage)){
            # yInner <- a[[1]]
            leverage[!(a[[1]] < maxBase & ( (a[[1]] %% roundBase) > 0))] <- FALSE
            if(any(preRounded)){
              if(any( preRounded & leverage)){
                leverage[preRounded] <- FALSE
                if(!any(leverage)){
                  warning("Due to preRounded, small inner cells may be revealed (easyCheck/leverageCheck)")
                }
              }
            }
            if(any(leverage)) {
              if(printInc) {
                cat(paste0("(Check:",sum(leverage),")")) 
                flush.console()
              }
              supRowsForce[leverage] = TRUE
              return_a <- FALSE
            } 
          }
        } 
      }
      
      if (return_a) {
        return(a)
      }
      
  }
  stop("Iteration limit exceeded")
}



PlsRoundSparseSingle  <- function(x,roundBase=3, yInner, yPublish = Matrix::crossprod(x,yInner)[,1,drop=TRUE],
                                     singleRandom = FALSE,
                                     suppPublish = yPublish < roundBase & yPublish > 0, #  Publiserte celler som skal undertrykkes (men må bruke iterasjon)
                                     yInnerExact = yInner,
                                     yPublishExact = yPublish,
                                  maxIterRows = 1000, supRowsForce = rep(FALSE, length(yInner)), identifyNew = TRUE, 
                                  preRounded = FALSE,
                                  step = 0, printInc = TRUE, plsWeights = NULL) {
  Pls1RoundHere <- get0("Pls1RoundFromUser", ifnotfound = Pls1Round) # Hack som gjør det mulig å bytte ut Pls1Round med annen algoritme

  printIncInput <- printInc
  
  roundBase <- as.integer(roundBase)


  if(identifyNew){
    suppInput <- yInner < roundBase & yInner > 0  #  Indre celler med verdier som er 'undertrykkbare'
    supRows <- (Matrix::rowSums(x[, suppPublish, drop = FALSE]) > 0 & suppInput) | supRowsForce 
    if(any(preRounded)){
      supRows[preRounded] <- FALSE
    }
  } else {
    supRows <- supRowsForce 
  }
  
  
  # To distribute correction when extreme correction caused by special input (preRounded)
  if (any(preRounded) &  any(supRows)) { #if (any(supRows)) {
    nRallOrdinary <- sum(yInner[supRows]%%roundBase)/roundBase
    nRallOrdinaryNot <- sum(supRows) - nRallOrdinary
    nRallOrdinaryCorrection <- (sum(yInnerExact) - sum(yInner))/roundBase
    correctionNotPart <- nRallOrdinaryCorrection/nRallOrdinaryNot
  } else {
    correctionNotPart <- 0
  }
  
  if (!singleRandom) if (sum(supRows) > maxIterRows) {
    randInd <- sample.int(sum(supRows), maxIterRows)
    supInds <- which(supRows)
    supRows[supRows] <- FALSE
    supRows[supInds[randInd]] <- TRUE
    printInc <- FALSE
    
    if(printIncInput){
      cat("#")
      flush.console()
    }
  }
  


  # Reduserer til antall rader som trengs
  bSupA <- x[supRows, , drop = FALSE]
  ySupp <- yInner[supRows]

  # Reduserer mer ved å fjerne unødvendike kolonner fungerer ikke på sparse ---- dessuten tar svært lang tid
  # cols1 =!duplicated(bSup,MARGIN=2)
  # bSupB = bSupA[,cols1,drop=FALSE]

  ## cols2 <- (Matrix::colSums(bSupA) > 0) & (Matrix::colSums(!bSupA) > 0)  # raskere beregning?

  colSumsbSupA <- Matrix::colSums(bSupA)
  cols2 <- (colSumsbSupA > 0) & ((NROW(bSupA)-colSumsbSupA)   > 0)  # raskere nå og mindre minnebruk.


  bSup <- bSupA[, cols2, drop = FALSE]

  yPublishCorrection <- yPublishExact[cols2] - yPublish[cols2]
  yPls <- t(as.matrix(Matrix::crossprod(bSup, Matrix(ySupp, ncol = 1))))
  
  plsWeightsHere <- plsWeights[cols2] 
  
  if (length(yPls) == 0 )   ## When 0 col in bSup
    singleRandom = TRUE
  
  if (sum(supRowsForce) & length(yPls) & length(ySupp)) {
    ySuppBase <- ySupp
    if(roundBase > 1L) {
      ySupp <- ySupp%%roundBase
    }
    ySuppBase <- ySuppBase - ySupp
    yPlsBase <- t(as.matrix(Matrix::crossprod(bSup, Matrix(ySuppBase, ncol = 1))))
    yPls <- yPls - yPlsBase
  } else {
    ySuppBase <- 0
  }
  
  
  #  yPls2 <- t(as.matrix(Matrix::crossprod(bSup, Matrix(ySupp, ncol = 1))))   # Not used
  
  
  correction <- TRUE  # -- For testing
  
  nR <- sum(ySupp)/roundBase
  
  if (correction) {
    yPls <- yPls + yPublishCorrection
    
    # Use old method when small difference
    if(correctionNotPart != 0)
      if(abs( (sum(yInnerExact) - sum(yInner))/roundBase) < 0.5001)
        correctionNotPart <- 0
      
    if(correctionNotPart == 0){   # Cannot be removed. correctionNotPart  set 0 also several lines above.
      nR <- round((sum(ySupp) + sum(yInnerExact) - sum(yInner))/roundBase) # old
    } else {
      nR <- round(nR + correctionNotPart*(length(ySupp)-nR))  # new  
    }
    
    nR <- min(nR, length(ySupp))
    nR <- max(nR, 0)
    
  } else {
    nR <- round(nR)
  }
  
  

  if (nR == 0 | singleRandom) {
    yR <- ySupp * 0L
    if (singleRandom){
      if(roundBase == 1L){
        makeFreq <- MakeFreq(matrix(sample.int(length(ySupp), nR, replace=TRUE), ncol=1))
        yR[makeFreq[,1]] <- makeFreq[,2]
      } else {
        yR[sample.int(length(ySupp), nR)] <- roundBase
      }
    }
  } else {
    yR <- Pls1RoundHere(bSup, ySupp, roundBase = roundBase, yPls = yPls, nR = nR, printInc=printInc, step = step, plsWeights=plsWeightsHere)
  }

  # Legger inn i ikke-reduserte data
  roundInner <- yInner
  roundInner[supRows] <- yR  + ySuppBase     ################################  ySuppBase   Ny
  roundPublish <- yPublish + Matrix::crossprod(bSupA, yR - ySupp)

  list(roundInner = roundInner, roundPublish = roundPublish, supRows = supRows)
}



Pls1Round <- function(x, y, roundBase = 3L, removeOneCols = FALSE, printInc = TRUE, yPls = NULL, nR = NULL, random = TRUE, wD=TRUE, step = 0, plsWeights) {
  # dgT med eller uten tApp/wD er muligheter ved for lite minne til roundBasecrossprod
  # wD fungerer raskt!!
  # cat("Pls1RoundFromUser")
  if(printInc) {cat("-"); flush.console()}
  if (is.matrix(x))
    x <- Matrix(x)  # Sparse matrix
  if (removeOneCols)
    x <- x[, (colSums(x) > 1), drop = FALSE]
  if (is.null(yPls))
    yPls <- t(as.matrix(Matrix::crossprod(x, Matrix(y, ncol = 1))))
  yR <- rep(0L, length(y))
  if (random)
    ind <- as.list(sample.int(length(y)))
  else ind <- as.list(seq_len(length(y)))
  
  indInv = vector("list",0)
  
  if (is.null(nR))
    nR <- round(sum(y)/roundBase)
  
  if(nR==0)
    return(yR)
  if(nR==length(y))
    return(rep(roundBase , length(y)))
  
  if(!is.null(plsWeights)){
    if(printInc) {cat("w"); flush.console()}
    plsWeights <- Diagonal(x=plsWeights)
    x <- x %*% plsWeights
    yPls <- yPls %*% plsWeights 
  }
  
  if(printInc) {cat("*"); flush.console()}
  #startTime <- Sys.time()
  #if(dgT){
  dgTBase <- As_TsparseMatrix(roundBase * Matrix::tcrossprod(x)) # dgTBase <- as(roundBase * Matrix::tcrossprod(x),"dgTMatrix") #flaskehals
    dgTi <- dgTBase@i +1L
    dgTj <- dgTBase@j +1L
    dgTx <- dgTBase@x
    rm(dgTBase)
    if(wD){
      dd = diff(dgTj)
      if(length(dd)>0){
        if(max(dd)>1 | min(dd)<0){
          warning("Not required sorting in dgTMatrix. Manual sorting will be done.")
          ord <- order(dgTj)
          dgTj <- dgTj[ord]
          dgTi <- dgTi[ord]
          dgTx <- dgTx[ord]
          dd = diff(dgTj)
        }
        wd <- c(1L,1L+which(dd==1L),length(dgTj)+1L)
      } else
        wd <- c(1L,2L)
      GetInd <- function(i,x){matlabColon(x[i],x[i+1L]-1L)}
    }
  #}
  #else Not used
  #  roundBasecrossprod <- as.matrix(roundBase * Matrix::tcrossprod(x))  # Much faster with as.matrix here
  #roundBasecrossprod <- as(roundBase * Matrix::tcrossprod(x),"dgTMatrix")  # Relativt treg
  if(printInc) {cat("*"); flush.console()}
  
  
  ind = as.integer(ind)
  indInv = integer(0)
  nBase = 0L
  
  coe <- Matrix::tcrossprod(x, yPls)
  
  if (roundBase != 1L) {
    UpBase <- function(returnIf1 = FALSE) {
      pf <- parent.frame()
      k <- which.max(coe[ind])
      ik <- ind[k]
      if (returnIf1) {
        if (k == 1) {
          return(FALSE)
        }
      }
      pf$yR[ik] <- roundBase
      pf$nBase <- nBase + 1L
      pf$indInv <- c(ind[k], indInv)
      pf$ind <- ind[-k]
      ii <- GetInd(ik, wd)
      ix <- dgTi[ii]
      pf$coe[ix] <- coe[ix] - dgTx[ii]
      TRUE
    }
    
    DownBase <- function(returnIf1 = FALSE) {
      pf <- parent.frame()
      k <- which.min(coe[as.integer(indInv)])
      ik <- indInv[k]
      if (returnIf1) {
        if (k == 1) {
          return(FALSE)
        }
      }
      pf$yR[ik] <- 0
      pf$nBase <- nBase - 1L
      pf$ind <- c(indInv[k], ind)
      pf$indInv <- indInv[-k]
      ii <- GetInd(ik, wd)
      ix <- dgTi[ii]
      pf$coe[ix] <- coe[ix] + dgTx[ii]
      TRUE
    }
  } else {   # roundBase == 1L
    UpBase <- function(returnIf1 = FALSE) {
      pf <- parent.frame()
      k <- which.max(coe[ind])
      ik <- ind[k]
      if (returnIf1) {
        if (k == 1) {
          return(FALSE)
        }
      }
      if (pf$yR[ik] == 0L) {
        pf$indInv <- c(ind[k], indInv)
      } else {
        pf$indInv <- c(ind[k], indInv[indInv != ind[k]])
      }
      ### pf$ind <- ind[-k]
      ii <- GetInd(ik, wd)
      ix <- dgTi[ii]
      pf$coe[ix] <- coe[ix] - dgTx[ii]
      pf$yR[ik] <- pf$yR[ik] + 1L
      pf$nBase <- nBase + 1L
      TRUE
    }
    
    DownBase <- function(returnIf1 = FALSE) {
      pf <- parent.frame()
      k <- which.min(coe[as.integer(indInv)])
      ik <- indInv[k]
      if (returnIf1) {
        if (k == 1) {
          return(FALSE)
        }
      }
      pf$yR[ik] <- pf$yR[ik] - 1L
      pf$nBase <- nBase - 1L
      ### pf$ind <- c(indInv[k],ind)
      if (pf$yR[ik] == 0L) {
        pf$indInv <- indInv[-k]
      }
      ii <- GetInd(ik, wd)
      ix <- dgTi[ii]
      pf$coe[ix] <- coe[ix] + dgTx[ii]
      TRUE
    }
  }
  
  
  step = max(1L, as.integer(step))
  
  i = 0L
  while(nBase<nR){
    i = i + 1L
    if (printInc)
      if (i%%max(1, round(nR/10)) == 0) {
        cat(".")
        flush.console()
      }
    UpBase()
    if(step > 1)
      if(i%%step == 0){
        doDown = TRUE
        for(i in seq_len(round(step/2))){
          if(doDown){ 
            #cat("D")
            doDown =DownBase(TRUE)  
          } 
        }
        #cat("|\n")
      }
  }
  for (i in 1:(nR+100)) {
    if (printInc)
      if (i%%max(1, round(nR/10)) == 0) {
        cat(":")
        flush.console()
      }
    doDown = TRUE
    for(i in seq_len(step)){
      if(doDown){ 
        #cat("d")
        doDown =DownBase(TRUE)   
      } 
    }
    #cat("|\n")
    if(nBase == nR){
      if(printInc) {cat("="); flush.console()}
      return(yR)
    }
    while(nBase<nR){
      UpBase()  
    }
  }
  if(printInc) {cat("="); flush.console()}
  yR
}















