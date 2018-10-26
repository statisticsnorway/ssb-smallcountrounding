#' function roundcube
#' 
#' 
#' This function rounds small counts in a set of hypercubes D produced by the function redcube 
#' and searches for a solution with smallest possible deviations from the original hypercube at
#' some aggregated levels.
#' 
#' @encoding UTF8
#'
#' @param rc The list of outpts from redcube
#' @param sort An ordered list of variables in hypercubes in D meant for priority sorting of the reduced 
#'             hypercube B before rounding. Not all variables in D should be included.
#' @param control A list of marginals of the hypercubes in D where deviations of aggregated rounded counts
#'                are checked against original counts.
#' @param minit Minimum number of searches to be carried out.
#' @param maxit Maximum number of searches to be carried out.
#' @param maxdiff If maximum difference in "control" is no larger than maxit, the stop search. 
#' @param seed Input seed for first systematic random search.
#'
#' @return 
#'    Ar: The rounded version of A
#'    Br: The rounded version of B
#'    D: The original hypercube of interest.
#'    Dr: The rounded version of D. The final table of interest.
#'    maxdiff: The largest absolute difference between cells D and Dr among cells in the control list.
#'    nmaxdiff: The number of occurences if Maxdiff
#'    
#' @keywords internal    
#' 
#' @importFrom stats runif
#' 
#' @author Johan Heldal, January 2018
#'
roundcube <- function(rc,sort,control,minit,maxit,maxdiff,seed) {
  A <- rc[[1]]  # The hypercube dataframe with all cells and variables..
  B <- rc[[2]]  # The reduced dataframe.
  C <- rc[[3]]  # The rows in A that are not in B.
  D <- rc[[4]]  # The tables/hypercubes of interest.
  Dr <-rc[[5]]  # The small counts (<b) in D.
  b <- rc[[6]]  # The rounding base.
  d <- rc[[7]]  # The list of variabel vectors defining D.
  nin <- rc[[8]]  # The count variable.
  nrB <- nrow(B)
  ncB <- ncol(B)
  Avars <- colnames(A[,1:(ncB-1)])
  #
  # Steplength for systematic sampling.
  #
  ss <- round(sum(B[,nin])/b)
  sumn <- sum(B[,nin])
  if (sumn>=b/2) {stp <- sumn/ss} 
  if (sumn<b/2)  {stp <- 0}
  #
  # Calculate control marginals for the reduced matrix B
  #
  S <- aggrtab(B,control,FALSE,nin=nin,nout=nin)[[2]]
  lS <- length(S)
  #
  # Initiate iterations
  #
  draw <- rep(0,nrB)
  set.seed(seed,kind="Mersenne-Twister")
  mmdiff <- 99999
  #
  # Start iterations
  #
  iter <- 0
  while ((iter <- iter+1) <= maxit) {
    #
    # For each iteration,sort the records in B by sort*u, u random. 
    # Result in Bu
    #
    u <- runif(nrB)
    Bu <- cbind(B,u)
    colnames(Bu) <- c(colnames(B),"u")
    sortvars <- c(sort,"u")        # Names of the sort variables
    nvsort <- length(sortvars)
    sortvars2 <- paste(rep("Bu$",nvsort),sortvars,sep="",collapse=",")
    txt <- paste("order(",sortvars2,")")
    Bu <- Bu[eval(parse(text=txt)),] # Sorting is carried out 
    #
    # Round: Draw a systematic pps-sample of cells in Bu with steplength stp.
    # Set the selected cell counts to b and the rest to 0.
    # Result in Bm.
    #
    csn <- cumsum(Bu[,ncB])
    strt <- runif(1)
    draw[1] <- strt*stp
    for (j in 1:(nrB-1)) {
      draw[j+1] <- draw[j] + stp*(csn[j]>= draw[j])
    }
    m <- b*(csn>=draw)         # Rounded reduced counts, this solution
    Bm <- cbind(Bu[,Avars],m)  # Rounded reduced cube, this solution
    #
    # Calculate the control marginals Sm for the rounded reduced cube
    #
    Sm <- aggrtab(Bm,control,FALSE,nin="m",nout="m")[[2]] # List of 
    Sk <- as.data.frame(NULL)  # "Current" original and rounded control table counts and differences. 
    Sd <- as.data.frame(NULL)  #   
    k <- 0 
#cat("iter = ",iter,"\n")   
    while((k <- k+1) <= lS) {
#cat("k = ",k,"\n")
      mm <- Sm[[k]]$m         # Frequencies of rounded reduced control table no. k.   
      Sk <- cbind(S[[k]],mm)  # Combine original and rounded counts in reduced control table no. k
      diff <- mm - S[[k]][,nin] # Differences between rounded and original controls no. k
      Sk <- cbind(Sk,diff)    
      ncSk <- ncol(Sk)
      vSk <- colnames(Sk)
      vSk <- c(vSk[(ncSk-2):ncSk],vSk[1:(ncSk-3)]) # Interchange columns
      Sk <- Sk[,vSk]                               # Put "nin", mm and diff first
      dk <- vSk[4:ncSk]
      if (k==1) {Sd <- Sk}
      if (k>1) {
        ncSd <- ncol(Sd)
        vSd <- colnames(Sd)
        dd <- vSd[4:ncSd]
        dd_k <- setdiff(dd,dk)  # colnames in dd not in dk
        dk_d <- setdiff(dk,dd)  # colnames in dk not in dd
        Sd_k <- as.data.frame(matrix(NA,nrow=nrow(Sd),ncol=length(dk_d)))
        colnames(Sd_k) <- dk_d
        Sd <- cbind(Sd,Sd_k)
        Sk_d <- as.data.frame(matrix(NA,nrow=nrow(Sk),ncol=length(dd_k)))
        colnames(Sk_d) <- dd_k
#cat("Sk_d = ","\n")
#print(Sk_d[1:20,])        
        Sk <- cbind(Sk,Sk_d)
        Sk <- Sk[,colnames(Sd)]
        Sd <- rbind(Sd,Sk)
        vSd <- colnames(Sd)
      } 
    }      
    mdiff <- max(abs(Sd$diff))    # Maximum diff this iteration
    ndiff <- sum(abs(Sd$diff)==mdiff) # No. of occurences of mdiff,this iteration
#cat("iter = ", iter, "mdiff =", mdiff, "ndiff = ", ndiff, "\n")    
    if (mdiff<mmdiff) {
      mmdiff <- mdiff              # Maximum difference of so far best solution 
      nmdiff <- ndiff             # No. of occurences of mmdiff in best soltion
      Br <- Bm                    # Br is best reduced solution "so far"
      Sr <- Sd                    # Control tables of the best solution
      cat("iter = ",iter," maxdiff = ",mmdiff," nmaxdiff = ",nmdiff,"\n")
    }
    if ((mdiff==mmdiff)&(ndiff<nmdiff)) {
      nmdiff <- ndiff
      Br <- Bm
      Sr <- Sd	  
      cat("iter = ",iter,", Maxdiff = ",mmdiff,", nmaxdiff = ",nmdiff,"\n")
    }
  }
# cat("colnames(C) = ", colnames(C), ", dim(C) = ", dim(C), "\n")
# cat("colnames(Br) = ", colnames(Br), ", dim(Br) = ", dim(Br),"\n")
  C <- C[,c(Avars,nin)]
  Br <- Br[,c(Avars,"m")]
  colnames(C) <- colnames(Br)
  C <- rbind(C,Br)
# cat("colnames(C) = ", colnames(C), ", dim(C) = ", dim(C),"\n","Avars = ", Avars, "\n")
  Cvars2 <- paste(rep("C$",length(Avars)),Avars,sep="",collapse=",")
  Ctxt <- paste("order(",Cvars2,")")
# cat("Cvars2 = ",Cvars2,"\n","Ctxt = ",Ctxt,"\n") 
  C <- C[eval(parse(text=Ctxt)),] # Sorting is carried out
  Avars2 <- paste(rep("A$",length(Avars)),Avars,sep="",collapse=", ")
  txt <- paste("order(",Avars2,")")
# cat("Avars2 = ",Avars2,"\n","txt = ",txt,"\n")   
  A <- A[eval(parse(text=txt)),] # Sorting is carried out
# cat("colnames(C) = ", colnames(C), ", dim(C) = ", dim(C),"\n")
# cat("colnames(Br) = ", colnames(Br), ", dim(Br) = ", dim(Br),"\n")
  m <- C$m
# cat("dim(A) = ", dim(A), ", length(m) = ", length(m),"\n")
  Ar <- cbind(A,m)
cat("colnames(Ar) = ", colnames(Ar), ", dim(Ar) = ", dim(Ar),"\n")
  Dr <- aggrtab(Ar,d,FALSE,nin="m",nout="m")[[2]]  
  return(list(Ar=Ar,Br=Br,D=D,Dr=Dr,Sm=Sm,Sr=Sr,maxdiff=mmdiff,nMaxdiff=nmdiff))  
}