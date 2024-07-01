

#' Function that returns a dataset 
#'
#' @encoding UTF8
#'
#' @param dataset Name of data set within the SmallCountRounding  package  
#' @param path When non-NULL the data set is read from "path/dataset.RData"
#'
#' @return The dataset
#' @export
#' @importFrom utils data
#' @importFrom SSBtools Hrc2DimList SSBtoolsData
#' 
#' @note Except for \code{"europe6"}, \code{"eHrc"}, \code{"eDimList"} and \code{"exPSD"}, the function returns the same datasets as \code{\link[SSBtools]{SSBtoolsData}}.
#' 
#' @seealso \code{\link[SSBtools]{SSBtoolsData}}, \code{\link[SSBtools]{Hrc2DimList}}
#'
#' @examples
#'  SmallCountData("z1")
#'  SmallCountData("e6")
#'  SmallCountData("eHrc")      #  TauArgus coded hierarchies
#'  SmallCountData("eDimList")  #  sdcTable coded hierarchies
#'  SmallCountData("exPSD")     #  Example data in presentation at Privacy in statistical databases
#'
SmallCountData <- function(dataset, path = NULL) {
  if (!is.null(path)) {
    filename <- paste(path, dataset, ".RData", sep = "")
    return(get(load(filename, envir = environment())))
  }
  if(dataset == "e6"){
    z <- data.frame(geo  = c("Iceland", "Portugal", "Spain"), 
                   eu = c("nonEU", "EU", "EU"),
                   year = rep(c("2018","2019"), each = 3),
                   freq = c(2,3,7,1,5,6), stringsAsFactors = FALSE)
    return(z)
  }
  if(dataset == "eHrc"){
    hierarchies = list(
      geo= c("EU", "@Portugal", "@Spain",  "nonEU",     "@Iceland"),
      year = c("2018", "2019"))
    return(hierarchies)
  }
  if(dataset == "eDimList"){
    return(Hrc2DimList(SmallCountData("eHrc"))) 
  }
  if(dataset == "exPSD"){
    return(data.frame(rows = rep(paste("row", 1:3, sep = ""), 5), 
                      cols = rep(paste("col", 1:5, sep = ""), each = 3), 
                      freq = c(6, 1, 0, 0, 2, 1, 1, 3, 1, 3, 1, 0, 4, 2, 2), 
                      stringsAsFactors = FALSE))
  }
  
  SSBtoolsData(dataset)
}
