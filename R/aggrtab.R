#' function aggrtab
#' 
#' @encoding UTF8
#' @noMd
#'
#' @param A A data frame representing a micro dataset or a frequency count hypercube. The (first) columns
#'          define the variables. If A is a hypercube the last column contains the number of units in each cell.
#' @param d A list d{[[j]]} whose elements are vectors of variable names from A defining marginal tables/cubes 
#'          of A that we are interested in.
#' @param micro Logical. TRUE if A is a micro dataset (default). FALSE if A i a frequency count hypercube.
#' @param nin Name of count variable if A is a hypercube. Default name: "n".
#' @param nout  Name of the frequency count variable in the output tables.
#'
#' @return
#'      D = {D[[j]]} of marginal tables/cubes of A spesified by the list d = d{[[j]]} generated 
#'      by aggregating over cells in A.
#'      
#' @keywords internal      
#' @importFrom  stats aggregate
#' 
#' @author Johan Heldal, November 2017
#'
aggrtab <- function(A,d,micro=TRUE,nin="n",nout="n") {
  #
  # Construct the vector Avars of all variables names entering the problem
  #
  d <- as.list(d)
  ld <- length(d)
  Avars <- d[[1]]
  k <- 1
  if (ld>1) {
    while ((k <- k+1) <= ld) {
      Avars <- union(Avars,d[[k]])
    }
  }
  ncA <- length(Avars)
  #
  # If A is a micro data set without cell counts, add a count column "n" with ones
  #
  if(micro==TRUE) {
    A <- A[,Avars]
    n <- rep(1,nrow(A))  # Add a column of ones as an aggregation variable
    A <- cbind(A,n)
    colnames(A) <- c(Avars,n)
    nin <- "n"
  }
  #
  # Aggregate to one cell per combination and set count variable name to "n"
  #
  A <- A[,c(Avars,nin)]
  A <- aggregate(A[,nin],by = as.list(A[,Avars]), FUN = "sum")
  colnames(A) <- c(Avars,nin)
  #
  # Sort A by the variables in Avars
  #
  Avars2 <- paste(rep("A$",ncA),Avars,sep="",collapse=", ")
  txt <- paste("order(",Avars2,")")
  A <- A[eval(parse(text=txt)),] # Sorting is carried out
  #
  # Generate the tables/cubes D defined by d.
  #
  D <- as.list(NULL)   # List of marginals of A to be built
  j <- 0
  while((j <- j + 1) <= ld) {
    dn <- d[[j]]               # Select the variable names for marginal table j
    dnn<- c(dn,nin)            # Add the cariable name "n" for the count variable
    D[[j]] <- A[,dnn]          # Selects the columns in A for marginal table j
    Dl <- as.list(as.data.frame(A[,dn]))
    #  
    # Aggregate A to marginal table j, D[[j]]
    #
    D[[j]] <- aggregate(D[[j]][,nin],by = Dl, FUN="sum")
    colnames(D[[j]]) <- c(dn,nout)
    #
    # Sort the cells in D[[j]] by the variables
    #
    ncd <- length(d[[j]])      # Number of variables in D[[j]]
    dn2 <- paste(rep("D[[j]]$",ncd),dn,sep="",collapse=", ")
    txt <- paste("order(",dn2,")")
    D[[j]] <- D[[j]][eval(parse(text=txt)),]    # Sorting is carried out
  }
  #
  return(list(Acube=A, Dcube=D, d=d, nin=nin, nout=nout))
}