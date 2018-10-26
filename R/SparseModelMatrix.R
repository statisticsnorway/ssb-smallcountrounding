#' Denne er erstattet av FormulaSums/ModelMatrix Overparameterized sparse model matrix with possible cross table
#'
#' All factor levels included. Constructed by calling fac2sparse() repeatedly instead of sparse.model.matrix().
#' Hierarchical variables handled when constructing cross table. Column names constructed from the cross table.
#'
#' @param formula A model formula
#' @param data data frame
#' @param mf model frame (alternative input instead of data)
#' @param hg List defining hierarchical variable groups
#' @param makeColnames Column names made when TRUE
#' @param crossTable Cross table in output when TRUE
#' @param total String used to name totals
#' @param sep String to separate when creating column names
#' @param printInc  Printing "..." to console when TRUE
#'
#' @return
#'   The model matrix or a list with model matrix and cross table.
#'
#' @keywords internal
#'
#' @examples
#'   z1 <- SmallCountData("z1")
#'   ModelMatrix(~region*hovedint,z1)
SparseModelMatrix = function(formula, data = NULL, mf = model.frame(formula, data = data), hg=HierarchicalGroups3(mf),
                             makeColnames=TRUE, crossTable=FALSE, total = "Total", sep="-", printInc=TRUE){
  hgid  = match(names(hg),colnames(mf))   # Merk hg endres ved endreing av input til HierarchicalGroups3

  hgcol = rep(0,NCOL(mf))
  for(i in seq_len(length(hg)))
    hgcol[hg[[i]]] = hgid[i]

  hgcoli = rep(0,NCOL(mf))
  for(i in seq_len(length(hg)))
    hgcoli[hg[[i]]] = i


  fac = attr(terms(as.formula(formula)),"factors")!=0
  faccol  = match(rownames(fac),colnames(mf))

  #print(colnames(fac)[k])


  #print(hgid)
  #print(faccol)
  #print(hgcol)


  firstROW = CharacterDataFrame(mf[1,hgid,drop=FALSE])

  firstROW = as.matrix(firstROW)

  firstROW[,] = total

  rownames(firstROW) = NULL

  allRows = firstROW

  m = fac2sparse(rep(1,NROW(mf)))

  nFac = NCOL(fac)

  for(k in seq_len(nFac)){
    if (printInc)
      if (k%%max(1, round(nFac/10)) == 0) {
        cat(".")
        flush.console()
      }
    ck = faccol[fac[,k]]
    #ur = UniqueRows(x[ ,ck,drop=FALSE]) # inn med fac2sparse(RowGroups... her


    if(makeColnames| crossTable){
      rg = RowGroups(mf[ ,ck,drop=FALSE],returnGroups  =TRUE)
      ur = rg[[2]]
      m = rbind(m,fac2sparse(rg[[1]])) # rBind(
      ur = CharacterDataFrame(ur)
      ur = as.matrix(ur)
      fr = firstROW[rep(1,NROW(ur)), ,drop=FALSE]
      #rownames(fr) = NULL
      fr[,hgcoli[ck]] = ur
      allRows = rbind(allRows,fr)
    } else
      m = rbind(m,fac2sparse(RowGroups(mf[ ,ck,drop=FALSE],returnGroups=FALSE))) # rBind(
  }
  if(makeColnames)
    rownames(m) <- MatrixPaste(allRows,sep=sep)

  if(!crossTable)
    return(Matrix::t(m))
  list(modelMatrix=m,crossTable=allRows)
}

#MatrixPaste = function(x, sep="_", forceCharacter=FALSE, stringEmpty = " "){
