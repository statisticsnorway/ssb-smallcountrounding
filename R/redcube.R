#' function redcube
#' 
#' This function produces a reduced small count frequency hypercube for rounding.
#' 
#' @encoding UTF8
#' @noMd
#' 
#' @param A A data frame representing a micro dataset or a frequency count hypercube. The (first) columns
#'          define the variables. If A is a hypercube the last column contains the number of units in each cell.
#'          If A is a micro dataset it is reduced to hypercube by the function aggrtab.
#' @param d A list d\{[[j]]\} whose elements are vectors of variable names from A defining marginal tables/cubes D
#'          of A that we are interested in.
#' @param b Rounding base. Counts in A less than b tat are contributing to counts less than b in the marginal 
#'          cubes D are selected from A. The selected dataframe is called B
#' @param micro Logical. TRUE if A is a micro dataset (default). FALSE if A i a frequency count hypercube.
#' @param nin Name of count variable if A is a hypercube. Default name: "n".
#'
#' @return
#'      A: The input dataframe reduced to a hypercube.
#'      
#'      B: The dataframe of small count rows selected from A.
#'      
#'      C: The dataframe of rows in A that are nor selected from A.
#'      
#'      D: The cubes defined by d.
#'      
#'      Dr: The small counts (<b) in D.
#'      
#'      The input elements d, b and nin
#'      
#' @keywords internal      
#' 
#' @author Johan Heldal, November 2017
#'
redcube <- function(A,d,b=3,micro=TRUE,nin="n") {
  AD <- aggrtab(A,d,micro=micro,nin=nin,nout=nin)
  A <- AD[[1]]
  ncA <- ncol(A) - 1       # Number of variables spanning A
  Avars <- colnames(A[,1:ncA])
  #
  #  For each j, create the tables Dr[[j]] of small cell counts in D[[j]].
  #  For each j identify the cells DA[[j]] in A aggregating to a cell in Dr[[j]] 
  #  Build the reduced hypercube B as the union of cells in {DA[[j]]}.
  #
  D <- AD[[2]]  
  ld <- length(d)
  Dr <- as.list(NULL)  # List of small count cells in D
  ninx <- paste(nin,".x",sep="",collapse=NULL)
  niny <- paste(nin,".y",sep="",collapse=NULL) 
  j <- 0
  while((j <- j + 1) <= ld) {
    dn <- d[[j]]               # Select the variable names for marginal table j
    Dr[[j]] <- D[[j]][D[[j]][,nin]<b,]              # Selects small cells in D[[j]]
    #
    # Select cells in A that aggregate to cells in Dr[[j]]
    #
    E <- merge(A,Dr[[j]],by.x=d[[j]],by.y=d[[j]])[,1:(ncA + 1)] 
    colnames(E) <- c(colnames(E[,1:ncA]),nin)
    
    if (j==1) {B <- E}
    if (j>1) {
      B <- merge(B,E,by=Avars,all=TRUE)
      B[,ninx][is.na(B[,ninx])] <- B[,niny][is.na(B[,ninx])]
    }
    B <- B[,1:(ncA+1)]
    colnames(B) <- c(colnames(B[,1:ncA]),nin) 
    #cat("j = ", j,"\n","dim(E) = ",dim(E),"\n","dim(D[[j]]) = ",dim(D[[j]]),"\n",
    #"dim(Dr[[j]]) = ",dim(Dr[[j]]),"\n","dim(B) = ",dim(B),"\n")
  }
  #
  # Calculate C, the dataframe og rows in A that are not in B
  #
  C <- merge(A,B,by=Avars,all.x=TRUE)
  C <- C[is.na(C[,niny]),]
  C <- C[,1:(ncA+1)]
  colnames(C) <- colnames(B)
  return(list(Acube=A, Bcube=B, Ccube=C, Dcubes=D, Drcubes=Dr, b=b, d=d, nin=nin))
}
