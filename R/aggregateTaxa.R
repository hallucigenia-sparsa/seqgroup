#' @title Aggregate a taxon abundance matrix
#'
#' @description Given a taxon and a lineage matrix, aggregate taxa belonging to the same
#' higher-level taxon at the specified taxonomic level by summing their abundance.
#'
#' @param x an abundance matrix with taxa as rows and samples as columns
#' @param lineages a lineage matrix where row names match row names of taxa; taxa must be provided in the same order as in x; first column gives kingdom, following columns give phylum, class, order, family, genus and species
#' @param taxon.level taxonomic level to which data set is aggregated (kingdom, phylum, class, order, family, genus or species)
#' @param unknown descriptors for unclassified taxa
#' @param keepUnknown keep taxa without higher-level taxon classification; cannot be used together with keepSumUnknown
#' @param keepSumUnknown keep the sum of taxa without higher-level taxon classification on the desired level as last row (row name: sumUnknown)
#' @param returnLineages return lineage matrix for higher-level taxon matrix; if true a list is returned with two elements: abundances and lineages
#' @return the aggregated abundance matrix; if returnLineages is true, a list with aggregated abundances and matching lineages
#' @examples
#' data("ibd_taxa")
#' data("ibd_lineages")
#' ibd_class=aggregateTaxa(x=ibd_taxa,lineages=ibd_lineages)
#' @export
aggregateTaxa<-function(x,lineages,taxon.level="class", unknown=c(NA,"unclassified"), keepUnknown = FALSE, keepSumUnknown=FALSE, returnLineages=FALSE){
  # this is needed, so lineages.higher.x can be filled without error caused by levels
  if(is.data.frame(lineages)){
    lineages=as.matrix(lineages)
  }
  if(nrow(lineages)!=nrow(x)){
    stop("The lineage matrix should have as many rows as the abundance matrix!")
  }
  if(keepSumUnknown && keepUnknown){
    stop("Please choose either keepUnknown or keepSumUnknown to keep unknown taxa.")
  }
  taxon.level=tolower(taxon.level) # conversion to lower case
  unknown.lineage=c("Dummy_phylum","Dummy_kingdom","Dummy_class","Dummy_order","Dummy_family","Dummy_genus","Dummy_species")
  unknown=c(unknown,unknown.lineage)
  indicesUnknown=c()
  levelIndex=NA
  if(taxon.level=="species"){
    levelIndex=7
  }else if(taxon.level=="genus"){
    levelIndex=6
  }else if(taxon.level=="family"){
    levelIndex=5
  }else if(taxon.level=="order"){
    levelIndex=4
  }else if(taxon.level=="class"){
    levelIndex=3
  }else if(taxon.level=="phylum"){
    levelIndex=2
  }else if(taxon.level=="kingdom"){
    levelIndex=1
  }else{
    stop("Requested taxon level not supported")
  }
  levelMembers=unique(lineages[,levelIndex])
  if(!keepSumUnknown && !keepUnknown){
    # discard unclassified or dummy higher-level taxa
    levelMembers=setdiff(levelMembers,unknown)
    rowNames=levelMembers
    rowNum=length(levelMembers)
  }else{
    # keep unclassified higher-level taxa, but do not count them as rows in the result matrix
    temp=setdiff(levelMembers,unknown)
    rowNames=temp
    rowNum=length(temp)
  }
  print(paste("Number of higher-level taxa:",length(levelMembers)))
  higher.x=matrix(NA,nrow=rowNum,ncol=ncol(x))
  lineages.higher.x=matrix(NA,nrow=rowNum,ncol=7)
  print(dim(lineages.higher.x))
  rownames(higher.x)=rowNames
  colnames(higher.x)=colnames(x)
  rowCounter=1
  for(levelMember in levelMembers){
    #print(levelMember)
    # get indices of member taxa
    member.indices=which(lineages[,levelIndex]==levelMember)
    if(length(member.indices)>0){
      #print(paste("level index:",levelIndex))
      #print(paste("row counter:",rowCounter))
      #print("first member index:")
      #print(member.indices[1])
      #print(lineages[member.indices[1],1:levelIndex])
      # pick the lineage of the first member, but only up to the selected level
      lineages.higher.x[rowCounter,1:levelIndex]=lineages[member.indices[1],1:levelIndex]
      if((keepUnknown || keepSumUnknown) && (levelMember %in% unknown)){
        indicesUnknown=c(indicesUnknown,member.indices)
      }else{
        # sum member taxa in each column
        vec=c()
        if(length(member.indices)==1){
          vec=colSums(t(as.matrix(x[member.indices,])))
          #print(dim(t(as.matrix(x[member.indices,]))))
        }else{
          vec=colSums(x[member.indices,])
          #print(dim(x[member.indices,]))
        }
        higher.x[rowCounter,]=vec
        rowCounter=rowCounter+1
      }
    # deal with missing values in taxonomic lineages
    }else{
      if(is.na(levelMember) && (keepSumUnknown || keepUnknown)){
        #print(paste("Keeping taxa with assignment NA for level",taxon.level))
        indicesUnknown=c(indicesUnknown,which(is.na(lineages[,levelIndex])))
      }else{
        warning(paste("No members found for higher-level taxon",levelMember))
      }
    }
  }
  indicesUnknown=unique(indicesUnknown)
  if(keepSumUnknown){
    print(paste("Adding sum of",length(indicesUnknown),"taxa with unknown classification."))
    unknownSum=colSums(x[indicesUnknown,])
    higher.x=rbind(higher.x,"sumUnknown"=unknownSum)
    lineages.higher.x=rbind(lineages.higher.x,"sumUnknown"=unknown.lineage)
  }else if(keepUnknown){
    print(paste("Appending",length(indicesUnknown),"taxa with unknown classification."))
    unknown.taxa.matrix=x[indicesUnknown,]
    higher.x=rbind(higher.x,unknown.taxa.matrix)
    lineages.higher.x=rbind(lineages.higher.x,lineages[indicesUnknown,])
  }
  if(returnLineages){
    rownames(lineages.higher.x)=rownames(higher.x)
    res=list(higher.x,lineages.higher.x)
    names(res)=c("abundances","lineages")
    return(res)
  }else{
    return(higher.x)
  }
}
