

### Fra RoundKostra

# This keeps integer but rowSums not
#RowSumsByApply = function(x){
#  apply(x,1,sum)
#}


#correctFormula <- function(x){
#  if(class(x)=="formula") return(x)
#  x <- trimws(as.character(x))
#  if(substr(x,1,1)!="~")
#    x = paste("~",x,sep="")
#  x
#}






#' Finding hierarchical variable groups
#'
#' As HierarchicalGroups() with eachName = TRUE, but output belonging to same mainName are combined.
#'
#' @param x Matrix or data frame containing the variables
#'
#' @return  List containing the groups.
#' @keywords internal
#' @export
#'
HierarchicalGroups2 <- function(x){
  a <- SSBtools::HierarchicalGroups(x,eachName = TRUE)
  b <- a[!duplicated(names(a))]
  for(i in 1:length(b)) b[[i]] = unique(unlist(a[names(a)==names(b)[i]]))
  b
}

#' Finding hierarchical variable groups
#'
#' As HierarchicalGroups() with eachName = FALSE, but output belonging to same mainName are combined.
#'
#' @param x Matrix or data frame containing the variables
#'
#' @return  List containing the groups.
#' @keywords internal
#' @export
#'
HierarchicalGroups3 <- function(x){
    a <- SSBtools::HierarchicalGroups(x,eachName = FALSE)
    b <- a[!duplicated(names(a))]
    for(i in 1:length(b)) b[[i]] = unique(unlist(a[names(a)==names(b)[i]]))
    b
  }


#' Make model formula from data taking into account hierarchical variables
#'
#' @encoding UTF8
#'
#' @param data data frame
#' @param hGroups Output from HierarchicalGroups2()
#' @param n Interaction level or 0 (all levels)
#' @param sim Include "~" when TRUE
#'
#' @return Formula as character string
#' @export
#'
#' @examples
#'  z3 <- SmallCountData("z3")
#'  MakeHierFormula(z3[,-7])
#'  MakeHierFormula(z3[,-7],n=2)
#'  MakeHierFormula(z3[,-7],n=0)
MakeHierFormula <- function(data=NULL,hGroups=HierarchicalGroups2(data),n=length(hGroups),sim=TRUE){
    if(n==0)
      sepS = ":"
    else
      sepS = "*"
    n = min(n,length(hGroups))
    m = AllNCombinations(sapply(hGroups,length),n)
    n = NROW(m)
    k = NCOL(m)
    x0 = rep("",k)
    z=rep("",n)
    for(i in seq_len(n)){
      mi = m[i,]
      x = x0
      for(j in seq_len(k))
        if(mi[j]) x[j] = hGroups[[j]][mi[j]]
      x = x[mi!=0]
      s = x[1]
      for(t in seq_len(length(x)-1))
        s = paste(s,x[t+1],sep=sepS)
      z[i] = s
    }
    if(!sim)
      return(paste(z,collapse=" + "))
    paste("~",paste(z,collapse=" + "),sep=" ")
  }



AllCombinations <- function(x = c(3,1,2),with0=TRUE,m=matrix(0,1,0)){
  if(!length(x)) return(m)
  nm <- NROW(m)
  AllCombinations(x[-1],with0,cbind(m[rep(seq_len(nm),x[1]+with0),],sort(rep(seq_len(x[1]+with0),nm))-with0))
}


AllNCombinations <- function(x = c(3,1,2),n=0,returnSorted=TRUE,returnList=FALSE){
  m = AllCombinations(x)
  rS = rowSums(m>0)
  if(n) return(m[rS==n, ,drop=FALSE])
  if(returnList){
    a <- vector("list",max(rS))
    for(i in seq_len(max(rS))) a[[i]] = m[rS==i, ,drop=FALSE]
    return(a)
  }
  m = m[!rS==0, ,drop=FALSE]
  rS = rS[!rS==0]
  if(returnSorted) return(m[order(rS), ,drop=FALSE])
  m
}


# Motsatt av easySdcTable:::MakeMicro
MakeFreq <- function(x,freqName="freq"){
  z = SSBtools::SortRows(x)
  uz = !duplicated(z)
  fr = matrix(-diff(c((NROW(x):1)[uz],0)),ncol=1,dimnames=list(NULL,freqName))
  b <- cbind(z[uz, ,drop=FALSE],fr)
  row.names(b) <- NULL
  b
}




### Fra Diversefunkjoner 


#' Combining columns of a matrix
#'
#' @param x  Matrix or vector
#' @param sep String used to combine columns
#' @param forceCharacter When FALSE single column input will keep to original class in output.
#' @param stringEmpty String used when input is empty (can be set to NULL)
#'
#' @return Character vector or possibly same vector as input
#' @details Each row in input will be combined to a single string using sep.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' MatrixPaste(matrix(1:12,3,4))
#' MatrixPaste(1:5)
#' MatrixPaste(1:5, forceCharacter=TRUE)
#' MatrixPaste(matrix(integer(0),3,0))
#' MatrixPaste(NULL)
#' }
MatrixPaste = function(x, sep="_", forceCharacter=FALSE, stringEmpty = " "){
  if(is.null(x)) return(stringEmpty)
  if(NCOL(x)==0) return(rep(stringEmpty,NROW(x)))
  if(NCOL(x)<=1){
    if(forceCharacter)
      return(as.character(as.vector(x)))
    return(as.vector(x))
  }
  apply(x , 1 , paste , collapse = sep )
}

#' @rdname MatrixPaste
#' @keywords internal
MatrixPaste1 = function(x,stringEmpty = "1") MatrixPaste(x,stringEmpty = stringEmpty)


