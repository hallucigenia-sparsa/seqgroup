#' @title Split samples into groups given metadata.
#' @description The function returns two groups situated in the low and high quantile of the given metadata item.
#' When two metadata items are provided, three splitting modi are available: congruent, complementary and inverse.
#' For groups with low and high values of a metadata item and two metadata items, congruent means that the low group
#' is low for both metadata items and the high group is high for both metadata items, inverse means that the low group is
#' low for the first metadata item and high for the second metadata item and vice versa for the high group whereas congruent means
#' that the first group is high for the first metadata item and the second high for the second metadata item.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param metadata a dataframe with metadata items as columns
#' @param metadata.name the name of a numeric metadata item to be used for sample splitting
#' @param metadata2.name the name of a second numeric metadata item to be used for sample splitting
#' @param mode the splitting mode; can be inverse, complementary or congruent; only relevant if a second metadata item is provided
#' @param quantile.def the thresholds on the lower and upper quantile which define the two sample groups (group 1 from 0 to first quantile, group 2 from second quantile to 1)
#' @return The function returns the abundances of group 1 and 2 (named group1 and group2) as well as the group-specific metadata values (named metadata1 and metadata2),
#' where second metadata item values can be empty.
#' @export
selectSamplesGivenMetadata<-function(abundances, metadata, metadata.name="", metadata2.name="", mode="congruent", quantile.def=c(0.1,0.9)){
  thresholds=quantile(metadata[[metadata.name]], quantile.def)
  metadata.values=metadata[[metadata.name]]
  metadata2.values=c()
  if(metadata2.name!=""){
    metadata2.values=metadata[[metadata2.name]]
  }
  indices.group1=which(metadata.values<thresholds[1]) # low quantile group
  indices.group2=which(metadata.values>thresholds[2]) # high quantile group
  if(metadata2.name!=""){
    thresholds2=quantile(metadata[[metadata2.name]], quantile.def)
    if(mode=="inverse"){
      indices.metadata2.group1=which(metadata2.values>thresholds2[2]) # group 1: high in metadata 2
      indices.metadata2.group2=which(metadata2.values<thresholds2[1]) # group 2: low in metadata 2
    }else if(mode=="complementary"){
      indices.metadata2.group2=which(metadata2.values>thresholds2[2])
    }else if(mode=="congruent"){
      indices.metadata2.group1=which(metadata2.values<thresholds2[1])
      indices.metadata2.group2=which(metadata2.values>thresholds2[2])
    }else{
      stop(paste("Mode",mode,"is not supported."))
    }
    if(mode=="complementary"){
      # metadata 1 should be high in first group
      indices.group1=indices.group2
      # metadata 2 should be high in second group
      indices.group2=indices.metadata2.group2
    }else{
      indices.group1=intersect(indices.group1,indices.metadata2.group1)
      indices.group2=intersect(indices.group2,indices.metadata2.group2)
    }
    if(length(indices.group1)==0){
      stop("No samples found in intersection of selected metadata for group 1.")
    }
    if(length(indices.group2)==0){
      stop("No samples found in intersection of selected metadata for group 2.")
    }
  }
  group1=abundances[,indices.group1]
  group2=abundances[,indices.group2]
  metadata.group1=metadata[indices.group1,]
  metadata.group2=metadata[indices.group2,]
  res=list(group1,group2, metadata.group1, metadata.group2)
  names(res)=c("group1","group2","metadata1","metadata2")
  return(res)
}

# Helper function to match age and gender for two data sets
# age1: age vector foor first data set
# gender1: gender vector for first data set
# age2: age vector for second data set
# gender2: gender vector for second data set
# range: allowed deviation for age
# TODO: to complete
matchAgeAndGender<-function(age1=c(), gender1=c(), age2=c(), gender2=c(), range=3){
  for(index in 1:length(age1)){
    queryage=age1[index]
    # try exact match first
    okindices=age2[age2==queryage]
    # if that fails, check for samples with age within the allowed range
    if(length(okindices)<1){
      okindices=c()
      for(secondindex in 1:length(age2)){
        if(age2[secondindex] < (queryage+range) && age2[secondindex] > (queryage-range)){
          okindices=c(okindices,secondindex)
        }
      }
    }
    print(paste("Found",length(okindices),"samples with age within range"))
    if(length(okindices>1)){
      # find matching gender
      #if(gender1)
    }
  }
}

