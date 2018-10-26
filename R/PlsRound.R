#' Small Count Rounding of Tabular Data
#'
#' Small count rounding via a dummy matrix and by an algorithm inspired by PLS
#'
#' Small count rounding of necessary inner cells are performed so that all small frequencies of cross-classifications to be published 
#' (publishable cells) are rounded. This is equivalent to changing micro data since frequencies of unique combinations are changed. 
#' Thus, additivity and consistency are guaranteed. 
#'
#' @encoding UTF8
#'
#' @param data Input data, data.frame or matrix
#' @param freqVar Variable holding counts (name or number)
#' @param formula Model formula defining publishable cells
#' @param roundBase roundBase
#' @param singleRandom Single random draw when TRUE (instead of algorithm)
#' @param crossTable	When TRUE, cross table in output and caculations via FormulaSums()
#' @param total	String used to name totals
#' @param maxIterRows See details
#' @param maxIter Maximum number of iterations
#' 
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
#' @seealso  See Round2() for rounding by other algorithm
#' @export
#'
#' @examples
#' RoundViaDummy(SmallCountData('z2'), 
#'               'ant', ~region + hovedint + fylke*hovedint + kostragr*hovedint, 10)
#' mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
#' a <- RoundViaDummy(SmallCountData('z3'), 'ant', mf, 5)
#' b <- RoundViaDummy(SmallCountData('sosialFiktiv'), 'ant', mf, 4)
#' print(cor(b[[2]]),digits=12) # Correlation between original and rounded
RoundViaDummy <- function(data, freqVar, formula, roundBase = 3, singleRandom = FALSE,
                          crossTable=TRUE, total = "Total",  maxIterRows = 1000, maxIter = 1E7) {
  cat("[")
  flush.console()
  previous_na_action <- options('na.action')
  options(na.action='na.pass')
  cat("{O")
  flush.console()
  if(crossTable){
    formulaSums <- FormulaSums(as.formula(formula), data = data, crossTable=TRUE,total=total,dropResponse=TRUE)
    x <- formulaSums$modelMatrix
    crossTab <- formulaSums$crossTable
    formulaSums <- NULL
  }
  else
    x <- ModelMatrix(as.formula(formula), data = data, sparse = TRUE)
  cat("}")
  flush.console()
  options(na.action=previous_na_action$na.action)

  #startTime <- Sys.time()
  if(anyNA(x))      # With na.action='na.pass' and sparse = TRUE no NA's will be produced now (zeros instead) - package ‘Matrix’ version 1.2-11
    x[is.na(x)] =0  # The code here is made for possible change in later versions
  #print(difftime(Sys.time(),startTime))

  yInner <- data[, freqVar]


  yPublish <- Matrix::crossprod(x, yInner)[, 1, drop = TRUE]
  a <- PlsRoundSparse(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish, singleRandom = singleRandom,maxIter=maxIter, maxIterRows=maxIterRows)
  cat("]\n")
  flush.console()


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


#' Overparameterized model matrix
#'
#' All factor levels included
#'
#' @param formula formula
#' @param data data frame
#' @param mf model frame (alternative input instead of data)
#' @param allFactor When TRUE all variables are coerced to factor
#' @param sparse When TRUE sparse matrix created by sparse.model.matrix()
#' @param formulaSums When TRUE, sparse matrix via FormulaSums()
#'
#' @return model matrix created via model.matrix(), sparse.model.matrix() or FormulaSums()
#' @importFrom stats model.frame model.matrix
#' @importFrom Matrix sparse.model.matrix
#' @export
#' @keywords internal
#'
#' @examples
#'   z1 <- SmallCountData("z1")
#'   ModelMatrix(~region*hovedint,z1)
ModelMatrix <- function(formula, data = NULL, mf = model.frame(formula, data = data), allFactor = TRUE, sparse = FALSE, formulaSums=FALSE) {
  if(formulaSums)
    return(FormulaSums(formula, data = data,
                                  makeNames=TRUE, crossTable=FALSE, total = "Total", printInc=TRUE,
                                  dropResponse = TRUE))
  for (i in 1:length(mf)) {
    if (allFactor)
      mf[[i]] <- as.factor(mf[[i]])
    if (is.factor(mf[[i]]))
      mf[[i]] <- AddEmptyLevel(mf[[i]])
  }
  if (sparse)
    return(sparse.model.matrix(formula, data = mf))
  model.matrix(formula, data = mf)
}


AddEmptyLevel <- function(x) factor(x, levels = c("tullnull", levels(x)))


#' PlsRoundSparse
#'
#' Avrunder basert på en algoritme inspirert av PLS-regresjon som forutsetter dummy matrise (Model matrix)
#'
#' @param x Model matrix
#' @param roundBase roundBase
#' @param yInner    inner cells
#' @param yPublish  cells to be published
#' @param singleRandom Single random draw when TRUE
#' @param yInnerExact Original yInner (when iteration)
#' @param yPublishExact Original yPublish (when iteration)
#'
#' @return rounded versions of yInner and yPublish
#' @importFrom  Matrix Matrix
#' @importFrom  methods as
#' 
#' @keywords internal
PlsRoundSparse <- function(x, roundBase = 3, yInner, yPublish = Matrix::crossprod(x, yInner)[, 1, drop = TRUE],
                           singleRandom = FALSE, maxIter = 1E6, maxIterRows = 1000) { # maxIter henger sammen med maxIterRows

  yInnerExact <- yInner
  yPublishExact <- yPublish


  i = 0
  while (i<maxIter) {
    i = i+1
    if (i == 1)
      a <- PlsRoundSparseSingle(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish,
                                singleRandom = singleRandom, yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows)
    else
      a <- PlsRoundSparseSingle(x = x, roundBase = roundBase, yInner = a[[1]], yPublish = a[[2]][, 1, drop = TRUE],
                                singleRandom = singleRandom, yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows)
      # suppRoundPublish = roundPublish<roundBase & roundPublish>0
      suppRoundPublish <- a[[2]] < roundBase & a[[2]] > 0
      if (!any(suppRoundPublish))
        return(a)
  }
  stop("Iteration limit exceeded")
}



PlsRoundSparseSingle  <- function(x,roundBase=3, yInner, yPublish = Matrix::crossprod(x,yInner)[,1,drop=TRUE],
                                     singleRandom = FALSE,
                                     suppPublish = yPublish < roundBase & yPublish > 0, #  Publiserte celler som skal undertrykkes (men må bruke iterasjon)
                                     yInnerExact = yInner,
                                     yPublishExact = yPublish,
                                  maxIterRows = 1000) {
  Pls1RoundHere <- get0("Pls1RoundFromUser", ifnotfound = Pls1Round) # Hack som gjør det mulig å bytte ut Pls1Round med annen algoritme

  roundBase = as.integer(roundBase)


  suppInput <- yInner < roundBase & yInner > 0  #  Indre celler med verdier som er 'undertrykkbare'

  supRows <- Matrix::rowSums(x[, suppPublish, drop = FALSE]) > 0 & suppInput

  printInc <- TRUE

  if(!singleRandom)
  if(sum(supRows)>maxIterRows){
    randInd = sample.int(sum(supRows),maxIterRows)
    supInds = which(supRows)
    supRows[supRows] = FALSE
    supRows[supInds[randInd]] = TRUE
    printInc <- FALSE
    {cat("#"); flush.console()}
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
  correction <- TRUE  # -- For testing
  if (correction) {
    yPls <- yPls + yPublishCorrection
    nR <- round((sum(ySupp) + sum(yInnerExact) - sum(yInner))/roundBase)
  } else nR <- round(sum(ySupp)/roundBase)

  
  if(length(yPls)==0)   ## When 0 col in bSup
    singleRandom = TRUE


  if (nR == 0 | singleRandom) {
    yR <- ySupp * 0L
    if (singleRandom)
      yR[sample.int(length(ySupp), nR)] <- roundBase
  } else yR <- Pls1RoundHere(bSup, ySupp, roundBase = roundBase, yPls = yPls, nR = nR, printInc=printInc)

  # Legger inn i ikke-reduserte data
  roundInner <- yInner
  roundInner[supRows] <- yR
  roundPublish <- yPublish + Matrix::crossprod(bSupA, yR - ySupp)

  list(roundInner = roundInner, roundPublish = roundPublish)
}





Pls1Round <- function(x, y, roundBase = 3L, removeOneCols = FALSE, printInc = TRUE, yPls = NULL, nR = NULL, random = TRUE,dgT=TRUE, wD=TRUE) {
  # dgT med eller uten tApp/wD er muligheter ved for lite minne til roundBasecrossprod
  # wD fungerer raskt!!
  #cat("Pls1RoundFromUser")
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

  if(printInc) {cat("*"); flush.console()}
  #startTime <- Sys.time()
  if(dgT){
    dgTBase <- as(roundBase * Matrix::tcrossprod(x),"dgTMatrix") #flaskehals
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
  }
  else
    roundBasecrossprod <- as.matrix(roundBase * Matrix::tcrossprod(x))  # Much faster with as.matrix here
  #roundBasecrossprod <- as(roundBase * Matrix::tcrossprod(x),"dgTMatrix")  # Relativt treg
  if(printInc) {cat("*"); flush.console()}
  for (i in 1:nR) {
    if (printInc)
      if (i%%max(1, round(nR/10)) == 0) {
        cat(".")
        flush.console()
      }
    if (i > 1){
      ii = GetInd(ik,wd)
      ix = dgTi[ii]
      coe[ix] <- coe[ix] - dgTx[ii]
    }
    else
      coe <- Matrix::tcrossprod(x, yPls)
    k <- which.max(coe[as.integer(ind)])
    ik <- ind[[k]]
    yR[ik] <- roundBase
    indInv <- c(ind[k],indInv)
    ind[k] <- NULL
  }
  absminmaxA = Inf
  absminmaxB = Inf
  for (i in 1:(nR+100)) {
    if (printInc)
      if (i%%max(1, round(nR/10)) == 0) {
        cat(":")
        flush.console()
      }
    ii = GetInd(ik,wd)
    ix = dgTi[ii]
    coe[ix] <- coe[ix] - dgTx[ii]
    k <- which.min(coe[as.integer(indInv)])
    if(k==1){
      if(printInc) {cat("="); flush.console()}
      return(yR)
    }

    ik <- indInv[[k]]
    yR[ik] <- 0
    ind <- c(indInv[k],ind)
    indInv[k] <- NULL
    ii = GetInd(ik,wd)
    ix = dgTi[ii]
    coe[ix] <- coe[ix] + dgTx[ii]
    k <- which.max(coe[as.integer(ind)])
    ik <- ind[[k]]
    yR[ik] <- roundBase
    indInv <- c(ind[k],indInv)
    ind[k] <- NULL
  }

  if(printInc) {cat("="); flush.console()}
  #endTime <- Sys.time()
  #print(difftime(endTime,startTime ))
  yR
}



















