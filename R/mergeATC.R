#' @title Merge ATC codes
#'
#' @description ATC codes are a classification system for therapeutic chemicals, see
#' https://en.wikipedia.org/wiki/Anatomical_Therapeutic_Chemical_Classification_System
#'
#' @param data matrix with only zeros or ones, with ATC codes as rows and samples as columns
#' @param level ATC level to which data should be merged
#' @return matrix with merged ATC codes
#' @export
mergeATC<-function(data, level=1){

  levels=c()
  # level 2: first 3 chars are used
  # level 3: first four chars
  # level 4: first five chars
  if(level>=2){
    level=level+1
  }
  # first seven chars are used
  if(level==5){
    level=level+2
  }

  # collect levels
  for(index in 1:nrow(data)){
    atccode=rownames(data)[index]
    code=paste(strsplit(atccode,split="")[[1]][1:level],collapse="")
    if(!(code %in% levels)){
      levels=c(levels,code)
    }
  }

  mergedatc=c()
  for(levelName in levels){
    print(paste("level",levelName))
    levelValues=c()
    # loop over samples (patients)
    for (colIndex in 1:ncol(data)){
      value=0
      # loop over rows
      for (rowIndex in 1:nrow(data)){
        atccode=rownames(data)[rowIndex]
        code=paste(strsplit(atccode,split="")[[1]][1:level],collapse="")
        # row belongs to current level
        if (code==levelName){
          if(data[rowIndex,colIndex]==1){
            value=1
            #print(paste(rownames(data)[rowIndex]," is part of level ",level,sep=""))
          }
        }
      } # end row loop
      levelValues=append(levelValues,value)
    } # end column loop
    mergedatc=rbind(mergedatc,levelValues)
  } # end level loop
  rownames(mergedatc)=levels
  colnames(mergedatc)=colnames(data)
  return(mergedatc)
}
