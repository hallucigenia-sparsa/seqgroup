#' @title Convert metadata into numeric form
#'
#' @description Binary metadata items are converted into binary numeric metadata items (with 0/1).
#' Categoric metadata with more than two categories can be binarized, such that
#' each category is represented by a separate binary metadata item. Numeric metadata are kept as is.
#' Dates, when specified, are converted into days since the reference date.
#' Note that constant metadata are removed.
#'
#' @param metadata a dataframe
#' @param yes the symbol used for the first value in a binary metadata item (e.g. "Y")
#' @param no the symbol used for the second value in a binary metadata item (e.g. "N")
#' @param na.threshold remove metadata with more than the maximum allowed percentage of missing values
#' @param to.skip metadata to skip from binarization
#' @param binarize convert categoric metadata items into as many binary metadata items as there are categories (if false, metadata with more than 2 categories are removed)
#' @param date.items names of metadata items that represent dates
#' @param format the format used for the date items (an example date fitting the default format is 26/1/80)
#' @param referenceDate reference date used for conversion of dates into numbers (the reference date format is always d/m/Y)
#' @return a purely numeric dataframe
#' @export
metadataToNumeric<-function(metadata, yes="Y", no="N", na.threshold=100, to.skip=c(), binarize=TRUE, date.items=c(),format="%d/%m/%y", referenceDate="1/1/1900"){
  binarizedMetadata=list()
  metadata.to.remove=c()
  skip.from.conversion=c()
  # loop metadata items
  for(name in names(metadata)){
    #print(paste("Processing",name))
    # remove leading or trailing white spaces
    for(j in 1:nrow(metadata)){
      metadata[[name]][j]=trimws(metadata[[name]][j])
    }
    # convert date into numeric
    if(name %in% date.items){
      metadata[[name]]=as.numeric(as.Date(metadata[[name]],format=format)-as.Date(referenceDate,format="%d/%m/%Y"))
    }else{
      if(is.factor(metadata[[name]])){
        levels=levels(metadata[[name]])
        levels=setdiff(levels,NA) # remove NA
        # process binary metadata
        if(length(levels)==2){
          skip.from.conversion=c(skip.from.conversion,name)
          numLevels=as.numeric(levels) # gives warnings (invalid factor level, NA generated) in case levels are non-numeric
          if(0 %in% numLevels && 1 %in% numLevels){
            print(paste("Binary metadata item",name,"is already numeric."))
            metadata[[name]]=as.character(metadata[[name]])
          }else{
            metadata[[name]]=as.character(metadata[[name]])
            yes.indices=which(metadata[[name]]==yes)
            no.indices=which(metadata[[name]]==no)
            # convert into a numeric metadata item, keep NA values
            replacement=rep(NA,length(metadata[[name]]))
            replacement[yes.indices]=1
            replacement[no.indices]=0
            metadata[[name]]=replacement
          }
        }else if(length(levels)>2){
          #print("Categoric factor")
          metadata.to.remove=c(metadata.to.remove,name)
          # binarize categoric metadata
          if(binarize && !(name %in% to.skip)){
            # there are less categories than samples
            if(length(levels)<nrow(metadata)){
              categoryNames=c()
              # initialize category-specific metadata items as absent
              # levels are NA-free
              for(level in levels){
                categoryName=paste(name,level,sep="-")
                categoryNames=c(categoryNames,categoryName)
                binarizedMetadata[[categoryName]]=rep(0,nrow(metadata))
              }
              # set presence values
              for(level in levels){
                for(j in 1:nrow(metadata)){
                  value=metadata[[name]][j]
                  if(!is.na(value)){
                    categoryName=paste(name,value,sep="-")
                    binarizedMetadata[[categoryName]][j]=1
                  }else{
                    # set all categories to NA
                    for(categoryName in categoryNames){
                      binarizedMetadata[[categoryName]][j]=NA
                    }
                  }
                }
              }
            }else{
              warning(paste("Cannot binarize categoric metadata item",name,"because it has as many categories as samples."))
            }
          }
        } # only 1 category - will be removed below
      } # end factor metadata
    } # not a date
  } # loop metadata
  # append binarized metadata
  if(binarize && length(binarizedMetadata)>0){
    for(append in names(binarizedMetadata)){
      #print(append)
      #print(length(binarizedMetadata[[append]]))
      metadata[[append]]=binarizedMetadata[[append]]
    }
  }
  metadata.filtered=filterMetadata(metadata,toFilter = metadata.to.remove,na.threshold = na.threshold)
  return(metadata.filtered)
}
