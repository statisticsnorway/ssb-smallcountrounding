# OldVersions 


RoundViaDummy_Version_0.3.0 <- function(data, freqVar, formula = NULL, roundBase = 3, singleRandom = FALSE,
                          crossTable=TRUE, total = "Total",  maxIterRows = 1000, maxIter = 1E7,
                          x = NULL,  hierarchies = NULL, Version = "tull", ...) {
  
  cat("[")
  flush.console()
  
  
  if (is.null(x) & is.null(formula) & is.null(hierarchies)) {
    freqVarName <- names(data[1, freqVar, drop = FALSE])
    hierarchies <- FindHierarchies(data[, !(names(data) %in% freqVarName)])
    # stop('formula, hierarchies or x needed')
  }
  
  
  if (!is.null(hierarchies) & !is.null(formula)) 
    stop("formula combined with hierarchies is not implemented")
  
  
  if (!is.null(hierarchies) & !is.null(x)) 
    warning("hierarchies ignored when x is supplied")
  
  if (!is.null(hierarchies) & is.null(x)) {
    x <- Hierarchies2ModelMatrix(data = data, hierarchies = hierarchies, crossTable = crossTable, total = total, ...)
    crossTable <- x$crossTable
    x <- x$modelMatrix
  }
  
  ## code below is as before hierarchies introduced 
  
  if(!is.null(x) & !is.null(formula))
    warning("formula ignored when x is supplied")
  
  if(is.null(x)){
    if(length(total)>1){
      total <- total[1]
      warning("Only first element of total is used when formula input.")
    }
    previous_na_action <- options('na.action')
    options(na.action='na.pass')
    cat("{O")
    flush.console()
    if(crossTable){
      formulaSums <- FormulaSums(formula = as.formula(formula), data = data, crossTable=TRUE,total=total,dropResponse=TRUE)
      x <- formulaSums$modelMatrix
      crossTab <- formulaSums$crossTable
      formulaSums <- NULL
    }
    else
      x <- ModelMatrix_Old_Version(as.formula(formula), data = data, sparse = TRUE)
    cat("}")
    flush.console()
    options(na.action=previous_na_action$na.action)
  } else {
    if(!is.logical(crossTable))
      crossTab <- crossTable
    crossTable <- TRUE
  }
  
  
  #startTime <- Sys.time()
  if(anyNA(x))      # With na.action='na.pass' and sparse = TRUE no NA's will be produced now (zeros instead) - package ‘Matrix’ version 1.2-11
    x[is.na(x)] =0  # The code here is made for possible change in later versions
  #print(difftime(Sys.time(),startTime))
  
  yInner <- data[, freqVar]
  
  
  yPublish <- Matrix::crossprod(x, yInner)[, 1, drop = TRUE]
  a <- PlsRoundSparse_Version_0.3.0(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish, singleRandom = singleRandom,maxIter=maxIter, maxIterRows=maxIterRows)
  cat("]\n")
  flush.console()
  
  
  if(crossTable)
    return(list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
                yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]),
                crossTable = crossTab))
  
  list(yInner = IntegerCbind(original = yInner, rounded = a[[1]]),
       yPublish = cbind(original = yPublish, rounded = a[[2]][, 1, drop = TRUE]))
}






PlsRoundSparse_Version_0.3.0 <- function(x, roundBase = 3, yInner, yPublish = Matrix::crossprod(x, yInner)[, 1, drop = TRUE],
                           singleRandom = FALSE, maxIter = 1E6, maxIterRows = 1000) { # maxIter henger sammen med maxIterRows
  
  yInnerExact <- yInner
  yPublishExact <- yPublish
  
  
  i = 0
  while (i<maxIter) {
    i = i+1
    if (i == 1)
      a <- PlsRoundSparseSingle_Version_0.3.0(x = x, roundBase = roundBase, yInner = yInner, yPublish = yPublish,
                                singleRandom = singleRandom, yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows)
    else
      a <- PlsRoundSparseSingle_Version_0.3.0(x = x, roundBase = roundBase, yInner = a[[1]], yPublish = a[[2]][, 1, drop = TRUE],
                                singleRandom = singleRandom, yInnerExact = yInnerExact, yPublishExact = yPublishExact, maxIterRows=maxIterRows)
    # suppRoundPublish = roundPublish<roundBase & roundPublish>0
    suppRoundPublish <- a[[2]] < roundBase & a[[2]] > 0
    if (!any(suppRoundPublish))
      return(a)
  }
  stop("Iteration limit exceeded")
}



PlsRoundSparseSingle_Version_0.3.0  <- function(x,roundBase=3, yInner, yPublish = Matrix::crossprod(x,yInner)[,1,drop=TRUE],
                                  singleRandom = FALSE,
                                  suppPublish = yPublish < roundBase & yPublish > 0, #  Publiserte celler som skal undertrykkes (men må bruke iterasjon)
                                  yInnerExact = yInner,
                                  yPublishExact = yPublish,
                                  maxIterRows = 1000) {
  Pls1RoundHere <- get0("Pls1RoundFromUser", ifnotfound = Pls1Round_Version_0.3.0) # Hack som gjør det mulig å bytte ut Pls1Round med annen algoritme
  
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


Pls1Round_Version_0.3.0 <- function(x, y, roundBase = 3L, removeOneCols = FALSE, printInc = TRUE, yPls = NULL, nR = NULL, random = TRUE,dgT=TRUE, wD=TRUE) {
  # dgT med eller uten tApp/wD er muligheter ved for lite minne til roundBasecrossprod
  # wD fungerer raskt!!
  #cat("Pls1Round_Version_0.3.0 ")
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
    dgTBase <- As_TsparseMatrix(roundBase * Matrix::tcrossprod(x))  #dgTBase <- as(roundBase * Matrix::tcrossprod(x),"dgTMatrix") #flaskehals
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




ModelMatrix_Old_Version <- function(formula, data = NULL, mf = model.frame(formula, data = data), allFactor = TRUE, sparse = FALSE, 
                        formulaSums=FALSE, printInc = FALSE) {
  if(formulaSums)
    return(formula = FormulaSums(formula, data = data,
                                 makeNames=TRUE, crossTable=FALSE, total = "Total", printInc=printInc,
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



