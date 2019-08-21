#' @title List core taxa across groups
#' @description List the indices of all taxa that occur at least min.occ times in each group.
#'
#' @param x an abundance matrix where rows are taxa and columns samples
#' @param groups group membership vector with as many entries as abundances has samples
#' @param min.occ minimal occurrence to count a taxon as present in each group
#' @return indices of core taxa
#' @examples
#' abundances=cbind(c(1,2,2),c(0,1,1),c(0,0,3))
#' rownames(abundances)=c("taxon1","taxon2","taxon3")
#' rownames(abundances)[getCoreTaxa(abundances,groups=c(1,2,3))]
#' @export
getCoreTaxa<-function(x,groups=c(),min.occ=1){
  if(length(groups)==0){
    stop("Please provide the group membership vector.")
  }
  core.taxa=c()
  unique.groups=unique(groups)
  for(i in 1:nrow(x)){
    present=TRUE
    taxon.vec=x[i,]
    for(group in unique.groups){
      indices=which(groups==group)
      occurrences=which(taxon.vec[indices]>0)
      if(length(occurrences)<min.occ){
        present=FALSE
      }
    } # end loop over groups
    if(present){
      core.taxa=c(core.taxa,i)
    }
  } # end loop over taxon rows
  return(core.taxa)
}
