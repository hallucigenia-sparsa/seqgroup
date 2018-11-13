#' @title Table union
#'
#' @description Merge two matrices row- or column-wise using their row
#' or column names (row names by default). Expects two matrices and appends
#' columns of the second table to the first such that row names in the
#' second table match row names in the first.
#' Rows with non-matching names are discarded by default, but can
#' be kept. In this case, their values are set to zero
#' in the table from which they are absent.
#'
#' @param table1 a matrix
#' @param table2 a matrix
#' @param keep.nonmatch keep non-matching rows or columns, with values missing in the other table filled with zeros
#' @param byRow compare tables row-wise
#' @return the merged table
#' @export
mergeTables<-function(table1, table2, keep.nonmatch=FALSE, byRow=TRUE){
  res=tableMatcher(table1 = table1, table2=table2, keep.nonmatch=keep.nonmatch, byRow = byRow)
  if(!byRow){
    # compensate for the initial transpose
    return(t(cbind(res$table1,res$table2)))
  }
  return(cbind(res$table1,res$table2))
}

#' @title Table intersection
#'
#' @description Discard rows or columns that are not present
#' in both tables (columns by default). The row or column
#' names are compared to determine whether they match.
#'
#' @param table1 a matrix
#' @param table2 a matrix
#' @param byRow compare tables row-wise
#' @param keepSumNonMatched keep the sum of non-matching rows as last row (row name: sumUnique)
#' @return a list with the two tables with matching sample names
#'
#' @export
intersectTables<-function(table1, table2, byRow=FALSE, keepSumNonMatched=FALSE){
  res=tableMatcher(table1 = table1, table2=table2, keep.nonmatch=FALSE, byRow = byRow, keepSumNonMatched = keepSumNonMatched)
  if(!byRow){
    res.t=list(t(res$table1),t(res$table2))
    names(res.t)=c("table1","table2")
    return(res.t)
  }
  return(res)
}

# internal tableMatcher function
tableMatcher<-function(table1, table2, keep.nonmatch = FALSE, keepSumNonMatched = FALSE, byRow=FALSE){
  # merge tables by their column names
  if(!byRow){
    table1=t(table1)
    table2=t(table2)
  }

  if(length(unique(rownames(table1)))!=nrow(table1)){
    stop("Table 1 contains non-unique names!")
  }
  if(length(unique(rownames(table2)))!=nrow(table2)){
    stop("Table 2 contains non-unique names!")
  }

  # sort both tables alphabetically
  table1=table1[sort(rownames(table1)),]
  table2=table2[sort(rownames(table2)),]

  rownames1=rownames(table1)
  rownames2=rownames(table2)
  common=intersect(rownames1,rownames2)
  print(paste("Number of matching row names:",length(common)))
  discard1=setdiff(rownames1,common) # row names only present in table 1
  print("Unique row names in table 1")
  print(discard1)
  discard2=setdiff(rownames2,common) # row names only present in table 2
  print("Unique row names in table 2")
  print(discard2)
  indicesDiscard1=c()
  indicesDiscard2=c()
  # find indices of row names to be discarded
  for(discard in discard1){
    indicesDiscard1=append(indicesDiscard1, which(rownames1==discard))
  }
  for(discard in discard2){
    indicesDiscard2=append(indicesDiscard2, which(rownames2==discard))
  }
  uniqueTable1=table1[indicesDiscard1,]
  uniqueTable2=table2[indicesDiscard2,]
  #print("Unique rows in table 1:")
  #print(rownames(uniqueTable1))
  #print("Unique rows in table 2:")
  #print(rownames(uniqueTable2))
  table1=table1[setdiff(1:nrow(table1),indicesDiscard1),]
  table2=table2[setdiff(1:nrow(table2),indicesDiscard2),]
  if(keepSumNonMatched){
    sumOfUniques1=colSums(uniqueTable1)
    sumOfUniques2=colSums(uniqueTable2)
    table1=rbind(table1,"sumUnique"=sumOfUniques1)
    table2=rbind(table2,"sumUnique"=sumOfUniques2)
  }
  if(keep.nonmatch==TRUE){
    # append rows only present in table 1 and set them to zero in table 2
    if(length(indicesDiscard1)>0){
      dummy=matrix(0,nrow=length(indicesDiscard1), ncol=ncol(table2))
      rownames(dummy)=rownames(uniqueTable1)
      print("Append rows to table 1:")
      print(rownames(dummy))
      table1=rbind(table1,uniqueTable1)
      table2=rbind(table2,dummy)
    }
    # append rows only present in table 2 and set them to zero in table 1
    if(length(indicesDiscard2)>0){
      dummy=matrix(0,nrow=length(indicesDiscard2), ncol=ncol(table1))
      rownames(dummy)=rownames(uniqueTable2)
      print("Append rows to table 2:")
      print(rownames(dummy))
      table2=rbind(table2,uniqueTable2)
      table1=rbind(table1,dummy)
    }
  }
  res=list(table1,table2)
  names(res)=c("table1","table2")
  return(res)
}
