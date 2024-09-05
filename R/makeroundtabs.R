#' function makeroundtabs
#'
#' This function creates a set of consistently rounded frequency count tables or
#' hypercubes by means of a version of small count rounding.
#' 
#' @encoding UTF8
#' @noMd
#' 
#' @seealso Dependencies:  \code{\link{aggrtab}}, \code{\link{redcube}}, \code{\link{roundcube}}
#'
#' @param A A data frame representing a micro dataset or a frequency count hypercube. The (first) columns
#'          define the variables. If A is a hypercube the last column contains the number of units in each cell.
#'          If A is a micro dataset it is reduced to hypercube by the function aggrtab.
#' @param b Rounding base. Counts in A less than b tat are contributing to counts less than b in the marginal 
#'          cubes D are selected from A. The selected dataframe is called B
#' @param d A list d\{[[j]]\} whose elements are vectors of variable names from A defining marginal tables/cubes D
#'          of A that we are interested in.
#' @param micro Logical. TRUE if A is a micro dataset (default). FALSE if A i a frequency count hypercube.
#' @param sort An ordered list of variables in hypercubes in D meant for priority sorting of the reduced 
#'             hypercube B before rounding. Not all variables in D should be included.
#' @param control A list of marginals of the hypercubes in D where deviations of aggregated rounded counts
#'                are checked against original counts.
#' @param nin Name of count variable if A is a hypercube. Default name: "n".
#' @param nout Name of the frequency count variable in the output tables.
#' @param minit  Minimum number of searches to be carried out.
#' @param maxit Maximum number of searches to be carried out.
#' @param maxdiff If maximum difference in "control" is no larger than maxit, the stop search. 
#' @param seed Input seed for first systematic random search.
#'
#' @return
#'    Ar: The rounded version of A
#'    
#'    Br: The rounded version of B
#'    
#'    D: The original hypercube of interest.
#'    
#'    Dr: The rounded version of D. The final table of interest.
#'    
#'    maxdiff: The largest absolute difference between cells D and Dr among cells in the control list.
#'    
#'    nmaxdiff: The number of occurences if Maxdiff
#' 
#' @export
#' @keywords internal
#' 
#' @author Johan Heldal, January 2018
#'
makeroundtabs <- function(A,b=3,d,micro="TRUE",sort,control,nin="n",nout="n",minit=3,maxit=3,
                          maxdiff=5,seed) {
  #
  # Step 1: Create a hypercube and a reduced hypercube
  #
  rc <- redcube(A,d,micro=micro,nin=nin) 
  #
  # Step 2: Round the reduced hypercube and establish rounded cube.
  #
cat("Er i makerountabs","\n")
cat("names(rc) = ", names(rc),"\n")
cat("rc$Acube[1:10,] = ","\n")  
print(rc$Acube[1,10,]) 
cat("rc$Bcube[1:10,] = ","\n")
print(rc$Bcube[1,10,])
cat("rc$nin = ",rc$nin,"\n")
  roundedcube <- roundcube(rc,sort,control,minit,maxit,maxdiff,seed)
  return(roundedcube)
  
}
