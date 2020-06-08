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

# Helper function to match age and gender for two data sets.
# Age is matched first, and gender is matched in case there is
# more than one sample that matches age within given range.
# If range is smaller than 1, gender is not matched.
# age1: age vector for query data set
# gender1: optional gender vector for query data set
# age2: age vector for target data set
# gender2: optional gender vector for target data set; needed if gender1 is given
# range: allowed deviation for age in years
# The method returns the indices of the selected target samples.
matchAgeAndGender<-function(age1=c(), gender1=c(), age2=c(), gender2=c(), range=1){
  selected.target.indices=c()
  if(length(gender1)>1 && length(gender2)==0){
    stop("If you provide a query gender vector, please provide a target gender vector.")
  }
  if(range<=0 && length(gender1)>0){
    gender1=c()
    warning("In order to match age and gender, please provide a range larger 0.")
  }
  # loop query age vector
  for(query.index in 1:length(age1)){
    queryage=age1[query.index]
    # try exact match first
    okindices=age2[age2==queryage]
    newIndexFound=FALSE
    if(length(okindices)>0){
      for(okindex in okindices){
        # select only one match
        if(!(okindex %in% selected.target.indices) && !newIndexFound){
          newIndexFound=TRUE
          # no gender match required
          if(length(gender1)==0){
            selected.target.indices=c(selected.target.indices,okindex)
          }
        }
      } # end loop indices found
    } # end test indices found
    # if exact age match fails or if gender is provided or if all target indices were already selected, check for samples with age within the allowed range
    if(length(gender1)>0 || !newIndexFound){
      # check within range
      if(range>0){
        okindices=c()
        # collect target samples with age within range
        for(target.index in 1:length(age2)){
          if(age2[target.index] <= (queryage+range) && age2[target.index] >= (queryage-range)){
            okindices=c(okindices,target.index)
          }
        }
        print(paste("Found",length(okindices),"samples with age within range"))
        # no matching age found within range
        if(length(okindices)<1){
          warning(paste("No matching age found for sample",query.index," and range ",range,". Consider expanding the range."))
        }
        # find matching gender
        else if(length(gender1)>0){
          newIndexFoundWithGender=FALSE
          for(okindex in okindices){
            if(gender1[query.index]==gender2[okindex]){
              if(!(okindex %in% selected.target.indices) && !newIndexFoundWithGender){
                selected.target.indices=c(selected.target.indices,okindex)
                newIndexFoundWithGender=TRUE
              }
            }
          } # end loop indices
          if(!newIndexFoundWithGender){
            warning(paste("Did not find target sample with matching gender not found before for sample",query.index))
            for(okindex in okindices){
              if(!(okindex %in% selected.target.indices) && !newIndexFound){
                selected.target.indices=c(selected.target.indices,okindex)
                newIndexFound=TRUE
              }
            } # end loop indices
          }
          if(!newIndexFound){
            warning(paste("Did not find target sample for query sample",query.index, "with matching age not found before"))
          }
        # no need to match gender
        }else if(length(gender1)==0){
          for(okindex in okindices){
            if(!(okindex %in% selected.target.indices) && !newIndexFound){
              selected.target.indices=c(selected.target.indices,okindex)
              newIndexFound=TRUE
            }
          } # end loop indices
          if(!newIndexFound){
            warning(paste("Did not find target sample for query sample",query.index, "with matching age not found before"))
          }
        }
      } # end range larger 0; with 0 range, nothing else can be done
    } # end no new index found or gender matching enabled
    if(!newIndexFound){
      warning(paste("No matching age or no new matchig age found for sample",query.index," and range ",range,". Consider expanding the range."))
    }
  } # end loop over query indices
  return(selected.target.indices)
}

