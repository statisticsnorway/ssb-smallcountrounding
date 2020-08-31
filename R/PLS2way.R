#' Two-way table from PLSrounding output
#' 
#' Output from \code{\link{PLSrounding}} is presented as two-way table(s) in cases where this is possible.
#' A requirement is that the number of main dimensional variables is two. 
#' 
#' When parameter \code{"variable"} is \code{"code"}, output is coded as  \code{"#"} (publish), \code{"."} (inner) and \code{"&"} (both).
#'
#' @param obj Output object from \code{\link{PLSrounding}}
#' @param variable  One of \code{"rounded"} (default), \code{"original"}, \code{"difference"} or \code{"code"}.
#'
#' @return A data frame
#' 
#' @importFrom SSBtools Match Unstack
#' @export
#'
#' @examples
#' # Making tables from PLSrounding examples 
#' z <- SmallCountData("e6")
#' a <- PLSrounding(z, "freq", formula = ~eu * year + geo)
#' PLS2way(a, "original")
#' PLS2way(a, "difference")
#' PLS2way(a, "code")
#' PLS2way(PLSrounding(z, "freq", formula = ~eu * year + geo * year), "code")
#' eHrc2 <- list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019"))
#' PLS2way(PLSrounding(z, "freq", hierarchies = eHrc2))
PLS2way = function(obj,  variable = c("rounded", "original", "difference", "code")){
  
  variable <- match.arg(variable)
  
  if(class(obj) !=  "PLSrounded")
    stop("Input must be output from PLSrounding")
  
  if(NCOL(obj$inner) !=  5)
    stop("Number of main dimensional variables must be two")
  
  nInner = Nlevels(obj$inner[,1])*Nlevels(obj$inner[,2]) 
  
  if(Nlevels(obj$inner[,1:2]) != nInner){
    stop("Inner cells must correspond to a square table")
  }
  
  if(NROW(obj$inner) != nInner){
    stop("Duplicated inner cells")
  }
  
  matchPI = Match(obj$publish, obj$inner)
  
  both = rbind(obj$inner,obj$publish[is.na(matchPI), ])
  
  both$code = "#"
  both$code[seq_len(nInner)] = "."
  both$code[matchPI[!is.na(matchPI)]] = "&"
  
  nBoth = Nlevels(both[,1])*Nlevels(both[,2])
  if(Nlevels(both[,1:2]) != nBoth){
    stop("All cells (inner + publish) must correspond to a square table")
  }
  
  if(NROW(both) != nBoth){
    stop("Duplicated cells when looking at inner and publish")
  }
  
  levels1 = unique(c(sort(unique(obj$inner[,1])), sort(unique(both[,1]))))
  levels2 = unique(c(sort(unique(obj$inner[,2])), sort(unique(both[,2]))))
  
  ord = order(as.integer(factor(both[,1], levels = levels1)),as.integer(factor(both[,2], levels = levels2)))
  
  colnames_both = colnames(both)
  both = suppressWarnings(Unstack(both[ord, ],match(variable,colnames_both),2,blockVar=1, returnRowData = FALSE))
  
  cb1 = match(colnames_both[1], colnames(both))
  rownames(both) = both[, cb1]
  both = both[,-cb1]
  both
}



Nlevels = function(x){
  NROW(unique(x,MARGIN=1))
}