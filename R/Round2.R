#' Small count rounding by various methods
#'
#' Run \code{makeroundtabs()} or \code{RoundViaDummy2()}
#'
#'
#' @encoding UTF8
#'
#' @param data Input data.frame. Same as input parameter \code{A} to \code{makeroundtabs()}
#' @param control As input to \code{makeroundtabs()}
#' @param method One of "roundtabs" (default), "viadummy", "viadummyAll" (\code{allTerms=TRUE})
#' @param ... Other paramameters as input to \code{makeroundtabs()}
#'
#' @return Output from  \code{\link{makeroundtabs}} or \code{\link{RoundViaDummy2}}
#' @export
#' 
#' @keywords internal 
#'
#' @seealso \code{RoundKostra} (SSB internal), \code{\link{RoundViaDummy}}, \code{\link{Lists2formula}}, \code{\link{MakeControl}}, \code{\link{FindMaxDiff}}
#'
#' @examples
#' \dontrun{
#' z <- SmallCountData("sosialFiktiv")
#' d <- list(c("region","mnd") , c("hovedint","mnd2") , c("fylke","hovedint","mnd2") , 
#'           c("kostragr","hovedint","mnd"))
#' con <- MakeControl(d,z)
#' sor <- names(z)[c(4,5,3,2,1)]
#'
#' roundedA <-  Round2(data=z,b=3,d=d,micro=FALSE,sort=sor,control=con, nin="ant",nout="Rndtall",
#'                     minit=2,maxit=10,maxdiff=5,seed=123,method="roundtabs")
#' roundedB <-  Round2(data=z,b=3,d=d,micro=FALSE,sort=sor,control=con, nin="ant",nout="Rndtall",
#'                     minit=2,maxit=10,maxdiff=5,seed=123,method="viadummy")
#' #10 rows of rounded data
#' roundedA$Ar[1:10,]  #"roundtabs"
#' cbind(z,roundedB$yInner)[1:10,] #"viadummy"
#'
#' # recalculate maxdiff  nMaxdiff
#' dA <- FindMaxDiff(roundedA$Ar,con,"ant","m")
#' dB <- FindMaxDiff(z,con,roundedB$yInner[,1],roundedB$yInner[,2])
#'
#' # Formula from d and control
#' Lists2formula(d,con,z)
#'
#' # Formula from another d
#' d2 <-list(sor)
#' Lists2formula(d2,con,z)
#' Lists2formula(d2,con) # Without knowing data
#' Lists2formula(d2,data=z) # Without control
#' Lists2formula(d2,data=z) # Without control and data
#' }
Round2 <- function(data,control,...,method=c("roundtabs","viadummy","viadummyAll")){
  method <- match.arg(method)
  if(method=="roundtabs"){
    a <- makeroundtabs(A=data,control=control,...)
  }
  if(method=="viadummy"){
    a <- RoundViaDummy2dots(data=data,control=control,...)
  }
  if(method=="viadummyAll"){
    a <- RoundViaDummy2dots(data=data,control=control,allTerms=TRUE,...)
  }
  a
}



List2string <- function(x,sepWithin = ":", sepBetween = " + "){
  f1 <- function(x) paste(x,collapse = sepWithin)
  paste(sapply(x,f1),collapse = sepBetween)
}


AddSim <- function(x) paste("~",x,sep="")


#' Make suggested control - input to makeroundtabs()
#'
#' @param d as input to \code{makeroundtabs()}
#' @param data data.frame
#' @param level Interaction level
#'
#' @return list 
#' @export
#' @keywords internal
#'
MakeControl <- function(d,data,level=2){
  n = length(d)
  hg = vector("list", n)
  for(i in 1:n)
    hg[[i]] = HierarchicalGroups2(data[,d[[i]],drop=FALSE])
    for(i in 1:n){
      if(i==1)
        f <- as.formula(MakeHierFormula(hGroups=hg[[i]],n=level))
      else
        f <- update(f,paste("~ . + ",MakeHierFormula(hGroups=hg[[i]],n=level,sim=FALSE)))
    }
  strsplit(trimws(strsplit(as.character(f)[[2]],split="\\+")[[1]]),":")
}


#' Make formula from input parameters to makeroundtabs()
#'
#' @param d d
#' @param control control
#' @param data data
#'
#' @return formula as string
#' @export
#' @importFrom SSBtools HierarchicalGroups2 MakeHierFormula
#' @importFrom  stats update
#' @keywords internal
#'
Lists2formula <- function(d,control=NULL,data=NULL){
  if(is.null(data)){
    if(is.null(control))
      return(AddSim(List2string(d,sepWithin = "*")))
    f1 <- List2string(control,sepWithin = ":")
    f2 <- List2string(d,sepWithin = ":")
    return(AddSim(paste(f1,f2,sep=" + ")))
  }
  n = length(d)
  hg = vector("list", n)
  for(i in 1:n)
    hg[[i]] = HierarchicalGroups2(data[,d[[i]],drop=FALSE])

  #unique(lapply(hg,names))
  if(is.null(control)){
    for(i in 1:n){
      if(i==1)
        f <- MakeHierFormula(hGroups=hg[[i]])
      else
        f <- update(as.formula(f),paste("~ . + ",MakeHierFormula(hGroups=hg[[i]],sim=FALSE)))
    }
    return(AddSim(as.character(as.formula(f))[[2]]))
  }
  Lists2formula(unique(lapply(hg,names)),control=control)
}



# Variant som aksepterer overflødige parametere
RoundViaDummy2dots <- function(data,..., b, d, nin, micro=FALSE,control=NULL,allTerms=FALSE,
                               singleRandom = FALSE, useFormulaSums = FALSE)
  RoundViaDummy2(data=data,b=b,d=d,nin=nin,micro=micro,control=control,singleRandom=singleRandom)


#' RoundViaDummy2
#'
#' RoundViaDummy with input as makeroundtabs
#'
#' @param data data
#' @param b b
#' @param d d
#' @param nin nin
#' @param micro micro
#' @param control control
#' @param allTerms Use all interaction terms in formula instead of using control
#' @param singleRandom Single random draw when TRUE (instead of algorithm)
#'
#' @return Output from RoundViaDummy extended with  "yControl" "maxdiff"  "nMaxdiff" "formula"
#' @export
#' @keywords internal
#'
RoundViaDummy2 <- function(data, b, d, nin, micro=FALSE,control=NULL,allTerms=FALSE,singleRandom = FALSE){
  if(micro)
    data = MakeFreq(data,nin)
  if(allTerms)
      formula <- Lists2formula(d=d,control=NULL,data=data)
    else
      formula <- Lists2formula(d=d,control=control,data=data)
  rvd <- RoundViaDummy(data=data, freqVar=nin, formula=formula, roundBase=b, singleRandom=singleRandom)
  if(!is.null(control))
    rvd <- c(rvd,FindMaxDiff(data,control,rvd[[1]][,1],rvd[[1]][,2],datareturn=TRUE),formula=formula)
  else
    rvd <- c(rvd,FindMaxDiff1(rvd[[2]]),formula=formula)
  if(!micro)
    return(rvd)
  c(rvd,data=data)
}

#' Calculate maxdiff  nMaxdiff
#'
#' @param data data
#' @param control control
#' @param original original
#' @param rounded rounded
#' @param datareturn datareturn
#'
#' @return maxdiff and nMaxdiff
#' @export
#' @keywords internal
#'
FindMaxDiff <- function(data,control,original,rounded,datareturn=FALSE){
  #m = ModelMatrix(as.formula(AddSim(List2string(control))),data,formulaSums=TRUE)
  m = FormulaSums(formula = as.formula(AddSim(List2string(control))),data = data,makeNames=datareturn)
  if(length(original)==1)
    original <-  as.matrix(data[,original,drop=FALSE])
  if(length(rounded)==1)
    rounded <-  as.matrix(data[,rounded,drop=FALSE])
  yControl = as.matrix(Matrix::crossprod(m,cbind(original,rounded)))
  out = FindMaxDiff1(yControl)
  if(datareturn)
    return(c(list(yControl=yControl),out))
  out
}

FindMaxDiff1 <- function(x,doPrint=TRUE){
  x = cbind(x,difference=x[,2]-x[,1])
  absdiff = abs(x[,3])
  maxdiff = max(absdiff)
  indmax = which(absdiff==maxdiff)
  nmaxdiff = length(indmax)
  if(doPrint){
    cat("\n maxdiff = ",maxdiff," nmaxdiff = ",nmaxdiff,"\n\n")
  }
  list(maxdiff=maxdiff, nMaxdiff=nmaxdiff)
}

