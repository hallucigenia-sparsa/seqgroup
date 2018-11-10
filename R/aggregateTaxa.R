#' @title Aggregate a taxon abundance matrix
#'
#' @description Given a taxon and a lineage matrix, aggregate taxa belonging to the same
#' higher-level taxon at the specified taxonomic level by summing their abundance.
#'
#' @param x an abundance matrix with taxa as rows and samples as columns
#' @param lineages a lineage matrix where row names match row names of taxa; taxa must be provided in the same order as in x; first column gives kingdom, following columns give phylum, class, order, family, genus and species
#' @param taxon.level taxonomic level to which data set is aggregated (kingdom, phylum, class, order, family, genus or species)
#' @param unknown descriptors of unclassified taxa
#' @export
aggregateTaxa<-function(x,lineages,taxon.level="class", unknown=c(NA,"unclassified")){
  taxon.level=tolower(taxon.level) # conversion to lower case
  unknown=c(unknown,"Dummy_phylum","Dummy_kingdom","Dummy_class","Dummy_order","Dummy_family","Dummy_genus","Dummy_species")
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
  # discard unclassified or dummy higher-level taxa
  levelMembers=setdiff(levelMembers,unknown)
  #print(levelMembers)
  print(paste("Number of higher-level taxa:",length(levelMembers)))
  higher.x=matrix(NA,nrow=length(levelMembers),ncol=ncol(x))
  rownames(higher.x)=levelMembers
  colnames(higher.x)=colnames(x)
  rowCounter=1
  for(levelMember in levelMembers){
    # get indices of member taxa
    member.indices=which(lineages[,levelIndex]==levelMember)
    if(length(member.indices)>0){
      # sum member taxa in each column
      #print(levelMember)
      vec=c()
      if(length(member.indices)==1){
        vec=colSums(t(as.matrix(x[member.indices,])))
        #print(dim(t(as.matrix(x[member.indices,]))))
      }else{
        vec=colSums(x[member.indices,])
        #print(dim(x[member.indices,]))
      }
      higher.x[rowCounter,]=vec
    }else{
      warning(paste("No members found for",levelMember))
    }
    rowCounter=rowCounter+1
  }
  return(higher.x)
}
