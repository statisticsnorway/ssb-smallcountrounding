
ModelMatrixExtra = function(data, ..., formula = NULL, hierarchies = NULL, crossTable = FALSE, sparse=TRUE, viaOrdinary = FALSE, modelMatrix=NULL, dimVar=NULL){
  
  if(!is.null(modelMatrix)){
    
    sparseInput = class(modelMatrix)[1]!= "matrix"
    
    if(sparseInput & !sparse)
      modelMatrix = as.matrix(modelMatrix)
    
    if(!sparseInput & sparse)
      modelMatrix = Matrix(modelMatrix)
    
    if (!is.null(formula))
      warning("formula ignored when model matrix is supplied in input")
    
    if (!is.null(hierarchies)) 
      warning("hierarchies ignored when model matrix is supplied in input")
    
    if(is.logical(crossTable)){
      if(crossTable)
        warning('"crossTable=TRUE" ignored when model matrix is supplied in input. crossTable as data.frame input is possible.')
      return(modelMatrix)
    }
    return(list(modelMatrix=modelMatrix, crossTable=crossTable))
  }
  
  if(!is.null(dimVar) & is.null(hierarchies) & is.null(formula)){
    data = data[, dimVar]
  }
  
  if(viaOrdinary){
    previous_na_action <- options('na.action')
    options(na.action='na.pass')
    on.exit(options(na.action=previous_na_action$na.action))
    a <- ModelMatrix(data=data,..., formula = formula, hierarchies=hierarchies, crossTable =crossTable, sparse=sparse, viaOrdinary=viaOrdinary,modelMatrix=modelMatrix)
    if(is.list(a)){
      if(anyNA(a$modelMatrix))      
        a$modelMatrix[is.na(a$modelMatrix)] =0  
    } else {
      if(anyNA(a))      # With na.action='na.pass' and sparse = TRUE no NA's will be produced now (zeros instead) - package ‘Matrix’ version 1.2-11
        a[is.na(a)] =0  # Code here allow possible change in later versions
    }
  }
  
  ModelMatrix(data=data, ..., formula = formula, hierarchies = hierarchies, crossTable =crossTable, sparse=sparse, viaOrdinary=viaOrdinary,modelMatrix=modelMatrix)
}
