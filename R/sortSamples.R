#' @title Sort samples by group and optionally by time.
#'
#' @param abundances an abundance matrix with taxa as rows and samples as columns
#' @param groups group membership vector with as many entries as samples
#' @param time optional time vector with as many entries as samples; entries are numeric or characters of numbers
#' @param alphabetical option to order samples by group name in alphabetical order
#'
#' @return an index vector with the updated sample order
#' @export
sortSamples<-function(abundances, groups=c(), time=c(), alphabetical=FALSE){
  updated.index.order=c()
  unique.groups=unique(groups)
  if(alphabetical){
    unique.groups=sort(unique.groups)
  }
  # put samples in the same group next to each other
  for(group in unique.groups){
    group.member.indices=which(groups==group)
    if(length(time)>0){
      # collect time points for group members
      times.of.group.members=c()
      for(group.member.index in group.member.indices){
        times.of.group.members=c(times.of.group.members, as.numeric(as.character(time[group.member.index])))
      }
      #print("before sorting")
      #print(times.of.group.members)
      # sort in ascending order
      out=sort(times.of.group.members,index.return=TRUE)
      #print(times.of.group.members[out$ix])
      group.member.indices=group.member.indices[out$ix]
    } # time given
    updated.index.order=c(updated.index.order,group.member.indices)
  } # loop groups
  return(updated.index.order)
}

#' @title Select groups
#'
#' @param groups membership vector with as many entries as samples
#' @param selected.groups names of groups to keep
#' @param inverse keep groups that are not in selected groups
#' @param returnIndex return the indices of the selected samples
#' @return modified group membership vector
#' @export
selectGroups<-function(groups=c(), selected.groups=c(), inverse=FALSE, returnIndex=FALSE){
  updated.groups.vector=c()
  for(i in 1:length(groups)){
    group=groups[i]
    if((group %in% selected.groups && !inverse) || (!(group %in% selected.groups) && inverse)){
      if(returnIndex){
        updated.groups.vector=c(updated.groups.vector,i)
      }else{
        updated.groups.vector=c(updated.groups.vector,group)
      }
    }
  }
  return(updated.groups.vector)
}
