#' @title Filter metadata
#'
#' @description Remove metadata with more than the given percentage
#' of missing values, constant metadata and metadata with names as provided
#' in the toFilter vector.
#'
#' @param data a matrix where metadata are rows and samples are columns or a data frame of metadata items
#' @param toFilter names of metadata to remove
#' @param na.threshold percentage of NA values allowed, between 0 and 100
#' @return a data frame of filtered metadata items
#' @export
filterMetadata<-function(data, toFilter=c(), na.threshold=0){
  if(is.data.frame(data)){
    data=t(as.matrix(data))
  }
  onePerc=ncol(data)/100
  na.number.threshold=round(onePerc*na.threshold)
  print(paste("Allowed number of NAs:",na.number.threshold))
  indices.tokeep=c()
  for(i in 1:nrow(data)){
    values=unique(data[i,])
    values=setdiff(values,NA)
    numberNA=length(which(is.na(data[i,])))
    if(length(values)>1 && numberNA<=na.number.threshold){
      if(rownames(data)[i] %in% toFilter){
        print(paste("Skipping metadata item: ",rownames(data)[i],sep=""))
      }else{
        indices.tokeep=c(indices.tokeep,i)
      }
    }else{
      if(numberNA>na.number.threshold){
        print(paste("Filtering metadata",rownames(data)[i],"with index",i,"because it has more missing values (namely",numberNA,") than the allowed percentage"))
      }else{
        print(paste("Filtering metadata",rownames(data)[i],"with index",i,"because it is constant"))
      }
    }
  } # end loop
  return(as.data.frame(t(data[indices.tokeep,])))
}

#' @title Assign data types to metadata
#'
#' @description Treat all metadata with more than two values
#' as numeric, unless they are in the vector with categoric values.
#' Metadata with two values (without counting missing values) are treated as
#' binary (so categoric).
#' Metadata are supposed to have been filtered previously (so constant rows were removed).
#' Note that a message is printed to indicate which metadata type was assigned and to
#' warn about numeric metadata with non-numeric values.
#'
#' @param data a metadata matrix with rows as samples and columns as items or a dataframe
#' @param categoric metadata items to be treated as categoric data
#' @export
assignMetadataTypes<-function(data, categoric=c()){
  if(is.data.frame(data)==FALSE){
    data=as.data.frame(data)
  }
  for(i in colnames(data)){
    values=unique(data[[i]])
    values=setdiff(values,NA) # remove NA
    if(i %in% categoric){
      data[[i]]=factor(data[[i]],ordered=FALSE)
      print(paste("Metadata",i,"is categoric"))
    }else if(length(values)==2){
      data[[i]]=factor(data[[i]],ordered=FALSE)
      print(paste("Metadata",i,"is binary"))
    }else{
      print(paste("Metadata",i,"is numeric"))
      tryCatch(as.numeric(as.character(data[[i]])), warning=function(w) print(paste(i,"has non-numeric values!")))
      data[[i]]=as.numeric(as.character(data[[i]]))
    }
  } # loop column names
  return(data)
}

#' @title Replace missing values by group mean.
#'
#' @description Missing values are replaced either by the group mean
#' for numeric metadata or the most frequent group value for categoric
#' metadata. Missing values are only replaced if there are less
#' missing values than the threshold in the group, else the metadata item
#' concerned is removed. The output are metadata without missing values.
#'
#' @param metadata.df a data frame with metadata
#' @param groups a vector that specifies the group for each sample
#' @param na.threshold number of allowed missing values per group
#' @param metadata.to.skip metadata for which filling with the group mean does not make sense (e.g. date)
#' @export
setNAToGroupMean<-function(metadata.df, groups=c(), na.threshold=4, metadata.to.skip=c()){
  filter.indices=c()
  patientmean=NA
  # loop metadata items (which are columns)
  for(i in 1:length(names(metadata.df))){
    numberNA=length(which(is.na(metadata.df[[i]])))
    metadata.item=names(metadata.df)[i]
    if(numberNA==0){
      # nothing to do
    }else{
      if(numberNA<=na.threshold && length(which(metadata.to.skip==metadata.item))<1){
        # loop over samples (which are rows)
        for(j in 1:nrow(metadata.df)){
          if(is.na(metadata.df[j,i])){
            patient=groups[j]
            patient.indices=which(groups==patient)
            if(length(patient.indices)==1){
              # cannot replace NA, because only a single value is available
              print(paste("Metadata item",names(metadata.df)[i],"cannot be filled, because group",patient,"has only a single value!"))
              filter.indices=append(filter.indices,i)
            }else{
              # separate treatment for numeric and factor
              if(is.numeric(metadata.df[,i])){
                patientmean=mean(na.omit(metadata.df[patient.indices,i]))
              }else if(is.factor(metadata.df[,i])){
                values=unique(metadata.df[patient.indices,i])
                values=setdiff(values,c(NA))
                if(length(values)==1){
                  patientmean=values[1]
                  #print(paste("filled with",patientmean))
                }else{
                  print(paste("Metadata item",names(metadata.df)[i],"has more than one factor for group",patient,"!"))
                  patientmean=NA
                  filter.indices=append(filter.indices,i)
                }
              }
              metadata.df[j,i]=patientmean
              print(paste("Metadata item",names(metadata.df)[i],"is filled with group mean",patientmean,"for sample",rownames(metadata.df)[j]))
            }
          }
        } # loop over samples
      }# not too many missing values
      else{
        filter.indices=append(filter.indices,i)
        print(paste("Metadata item",names(metadata.df)[i],"has too many NA or is in the list of metadata to be skipped!"))
      }
    } # missing values present
  } # loop over metadata items
  keep.indices=setdiff(c(1:length(names(metadata.df))),filter.indices)
  return(metadata.df[,keep.indices])
}

# type can be either numeric, categoric, binary or catbin
# categoric does not include binary, but catbin includes both
getMetadataSubset<-function(metadata.df, type="numeric"){
  keep=c()
  for(i in 1:length(names(metadata.df))){
    if(type=="numeric" && is.numeric(metadata.df[,i])){
      keep=append(keep,i)
    }else if(is.factor(metadata.df[,i]) && (type=="categoric" || type=="binary" || type=="catbin")){
      if(type=="catbin"){
        keep=append(keep,i)
      }else{
        values=setdiff(unique(metadata.df[,i]),NA)
        if(type=="binary" && length(values)==2){
          keep=append(keep,i)
        }else if(type=="categoric" && length(values)>2){
          keep=append(keep,i)
        }
      }
    }
  }
  return(metadata.df[,keep])
}

# remove samples with indicated metadata item value
# if random is true, remove as many random samples as there are samples that have the particular metadata item value
# function returns filtered taxa and metadata objects in a list
filterSamplesWithMetadataValues<-function(taxa, metadata.df, metadata.item="", metadata.value=1, random=FALSE){
  metadata.values=metadata.df[[metadata.item]]
  to.filter=which(metadata.values==metadata.value)
  to.keep=setdiff(c(1:ncol(taxa)),to.filter)
  if(random){
    rand.indices=sample(c(1:ncol(taxa)))[1:length(to.keep)]
    to.keep=rand.indices
  }
  taxa=taxa[,to.keep]
  metadata.df=metadata.df[to.keep,]
  res=list(taxa, metadata.df)
  names(res)=c("taxa","metadata")
  return(res)
}

