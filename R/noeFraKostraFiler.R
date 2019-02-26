# Motsatt av easySdcTable:::MakeMicro
MakeFreq <- function(x,freqName="freq"){
  z = SSBtools::SortRows(x)
  uz = !duplicated(z)
  fr = matrix(-diff(c((NROW(x):1)[uz],0)),ncol=1,dimnames=list(NULL,freqName))
  b <- cbind(z[uz, ,drop=FALSE],fr)
  row.names(b) <- NULL
  b
}




### Fra Diversefunkjoner. Disse funksjonene legges også inn i SSBtools. Fjerner rd her


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


MatrixPaste1 = function(x,stringEmpty = "1") MatrixPaste(x,stringEmpty = stringEmpty)


