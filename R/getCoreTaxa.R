#' @title List core taxa across groups
#' @description List the indices of all taxa that occur at least min.occ times in each group.
#'
#' @param x an abundance matrix where rows are taxa and columns samples
#' @param groups group membership vector with as many entries as abundances has samples
#' @param min.occ minimal occurrence to count a taxon as present in each group
#' @param min.groupnum minimum occurrence across groups (fulfilling selected criteria within each group), by default in all
#' @param consecutive assumes samples in each group are ordered by time, if true, minimal occurrence has to be consecutive
#' @return indices of core taxa
#' @examples
#' data("floresgut_taxa")
#' data("floresgut_metadata")
#' groups=floresgut_metadata$host_subject_id
#' ct=getCoreTaxa(x=floresgut_taxa,groups=groups,consecutive=TRUE, min.occ=4)
#' rownames(floresgut_taxa)[ct]
#' @export
getCoreTaxa<-function(x,groups=c(),min.occ=1,min.groupnum=length(unique(groups)),consecutive=FALSE){
  if(length(groups)==0){
    stop("Please provide the group membership vector.")
  }
  core.taxa=c()
  unique.groups=unique(groups)
  groupnum=length(unique.groups)
  groupnumbers=rep(0,nrow(x)) # in the order of taxa
  # taxon loop
  for(i in 1:nrow(x)){
    present=TRUE
    taxon.vec=x[i,]
    # loop over groups
    for(group in unique.groups){
      indices=which(groups==group)
      occurrences=which(taxon.vec[indices]>0)
      if(length(occurrences)<min.occ){
        present=FALSE
      }
      # check whether occurrence is consecutive
      if(present && consecutive){
        searchStr=rep("n",length(indices))
        searchStr[occurrences]="y"
        # find longest consecutive stretch
        # exhaustive: 2 loops, avoid:
        # https://www.geeksforgeeks.org/maximum-consecutive-repeating-character-string/
        n=length(searchStr)
        currentCount=1
        count=0
        # traverse string
        for (j in 1:n){
          # if current char matches with next
          if(j < n && searchStr[j] == searchStr[j+1]){
            currentCount=currentCount+1
          }else{
            if(currentCount>count){
              count=currentCount
            }
            currentCount=1
          }
        } # go through search string
        # check if longest stretch is longer than minimum
        if(count<min.occ){
          present=FALSE
        }
      } # end check for consecutive occurrence
      # count presence across groups
      if(present){
        groupnumbers[i]=groupnumbers[i]+1
      }
    } # end loop over groups
    # present in all groups
    if(present){
      core.taxa=c(core.taxa,i)
    }
  } # end loop over taxon rows
  # taxa do not need to be present in all groups
  if(min.groupnum<groupnum){
    core.taxa=c()
    for(i in 1:length(groupnumbers)){
     if(groupnumbers[i]>=min.groupnum){
       core.taxa=c(core.taxa,i)
     }
    }
  }
  return(core.taxa)
}
