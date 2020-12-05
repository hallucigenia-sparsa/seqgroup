#' @title Bar plot of taxon composition with group support
#' @description Sort taxa by summed abundance across all samples and plot sorted taxon composition with a bar per sample
#' @details Note that taxa are always summed across all samples, also in the presence of a group membership vector, unless sumGroupwise is true.
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param groups group membership vector with as many entries as samples
#' @param aggregate if groups are given, plot the aggregate across the group (or its selected samples if randSampleNum is true); possible values: none, median and mean
#' @param taxon.color.map map of taxon-specific colors, should match row names; taxa not present in the color map will be colored in summedTaxonColor
#' @param group.color.map map of group-specific colors, should match group names
#' @param topTaxa number of top taxa to be plotted
#' @param sortGroupwise if true, samples are sorted according to groups
#' @param sumGroupwise if true, taxa are summed and sorted separately across samples within each group (if true, samples are always sorted group-wise)
#' @param group.order if a vector with group names (one for each group) is given, group samples will be sorted in the order indicated; can also be used to only plot selected groups
#' @param hide.taxa do not consider these taxa as top-abundant taxa, but keep them among Others
#' @param randSampleNum if larger 0, sortGroupwise is set to true and the indicated sample number is randomly selected for each group
#' @param summedTaxonColor the color of the summed taxa, by default gray
#' @param extendTaxonColorMap if true, taxa not in the taxon color map are added there and the extended color map is returned
#' @param legend add a legend with the color code
#' @param legend.shift increase/decrease this parameter to shift the color legend further to the right/left
#' @param legend.hidegroups do not show the group memberships in the legend
#' @param \\dots Additional arguments passed to plot()
#' @return if extendTaxonColorMap is true, the taxon color map is returned
#' @examples
#' data(ibd_taxa)
#' data(ibd_metadata)
#' groups=as.vector(ibd_metadata$Diagnosis)
#' # taxon abundances were prefiltered and therefore do not add up to 1
#' groupBarplot(ibd_taxa,groups=groups,randSampleNum=10)
#' # sum taxa group-wise for sorting instead across all samples
#' groupBarplot(ibd_taxa,groups=groups,sumGroupwise=TRUE, legend.hidegroups=TRUE)
#' data(ibd_lineages)
#' ibd_genera=aggregateTaxa(ibd_taxa,ibd_lineages,taxon.level = "genus")
#' groupBarplot(ibd_genera,groups=groups,randSampleNum=10)
#' @export
groupBarplot<- function(abundances, groups=c(), aggregate="none", taxon.color.map=NULL, group.color.map=NULL, topTaxa=10, sortGroupwise=TRUE, sumGroupwise=FALSE, group.order=c(), hide.taxa=c(), randSampleNum=NA, summedTaxonColor="#a9a9a9", extendTaxonColorMap=FALSE, legend=TRUE, legend.shift=1, legend.hidegroups=FALSE, ...){

  if(sumGroupwise && length(groups)==0){
    warning("Can only sum groupwise if groups are provided. Summing groupwise will be ignored.")
    sumGroupwise=FALSE
  }

  if(sumGroupwise){
    sortGroupwise=TRUE
  }

  if(sumGroupwise && length(hide.taxa)>0){
    stop("Summing group-wise cannot be used together with hiding taxa.")
  }

  if(length(groups)>0){
    if(length(groups)!=ncol(abundances)){
      stop("Each sample should have a group assigned.")
    }
  }

  if(aggregate != "none" && length(groups)==0){
    stop("For aggregation across groups, a group membership vector is needed.")
  }

  if(aggregate == TRUE){
    stop("Possible values for aggregate are none, median and mean.")
  }

  # if no groups are assigned, put all samples into the same default group
  if(length(groups)==0){
    groups=rep("all",ncol(abundances))
  }

  # assign taxon names if necessary
  if(is.null(rownames(abundances))){
    rownames(abundances)=as.character(1:nrow(abundances))
  }

  # assign sample names if necessary
  if(is.null(colnames(abundances))){
    colnames=as.character(1:ncol(abundances))
  }

  if(!is.na(randSampleNum) && randSampleNum>0){
    sortGroupwise=TRUE
  }else{
    randSampleNum=0
  }

  if(aggregate != "none"){
    sortGroupwise=TRUE
  }

  groupNum=length(unique(groups))
  prev.mar=par()$mar
  mar.scale=0.5 # maximal number of characters is scaled by this number to compute mar on the right side

  sorted=NULL
  superabundances=matrix(NA,nrow=1,ncol=1)
  taxon.colors=c()
  group.colors=c()

  group.indices.map=list()
  aggregated=matrix(NA,nrow=nrow(abundances),ncol=length(unique(groups)))
  rownames(aggregated)=rownames(abundances)
  counter=1

  # sort samples according to group membership
  if(sortGroupwise && (groupNum>1 || randSampleNum>0)){
    abundancesSorted=NULL
    #updated.groups=c()
    # loop groups
    for(group in unique(groups)){
      group.member.indices=which(groups==group)
      print(paste("Number of samples in group",group,":",length(group.member.indices)))
      if(randSampleNum>0){
        #print("Selecting random samples")
        if(length(group.member.indices)>randSampleNum){
          group.member.indices=sample(group.member.indices)[1:randSampleNum]
        }
      }
      group.indices.map[[group]]=group.member.indices
      if(!sumGroupwise){
        if(aggregate == "median"){
          group.indices.map[[group]]=counter
          aggregated[,counter]=apply(abundances[,group.member.indices],1,median) # compute median row-wise
        }else if(aggregate == "mean"){
          group.indices.map[[group]]=counter
          aggregated[,counter]=apply(abundances[,group.member.indices],1,mean) # compute mean row-wise
        }
      }
      counter = counter +1
    } # end loop groups
    indices=c()
    # for a single group or no group, order is irrelevant
    if(length(group.order)>1){
      for(group in group.order){
        indices=c(indices,group.indices.map[[group]])
      }
    }else{
      indices=unlist(group.indices.map)
    }
    if(aggregate == "none" || sumGroupwise){
      abundances=abundances[,indices]
      groups=groups[indices]
    }else{
      abundances=aggregated[,indices]
      if(length(group.order)>0){
        groups=group.order
      }else{
        groups=unique(groups)
      }
      colnames(abundances)=groups
    } # aggregate is not none
  } # sort group-wise

  hidden.taxa=NULL
  if(length(hide.taxa)>0){
    indices.hidden=c()
    for(hidden.taxon in hide.taxa){
      indices.hidden=c(indices.hidden,which(rownames(abundances)==hidden.taxon))
    }
    hidden.taxa=as.matrix(abundances[indices.hidden,])
    if(length(indices.hidden)==1){
      hidden.taxa = t(hidden.taxa)
    }
    keep.indices=setdiff(c(1:nrow(abundances)),indices.hidden)
    abundances=abundances[keep.indices,]
  }

  # sum taxa group-wise
  if(sumGroupwise && groupNum>1){
      # superabundances has each row as many times as there are groups; +1 for the summed taxon, which is group-specific
      # this allows to re-sort rows for each group separately but to keep the original colors (by duplicating colors)
      # values for samples not in the group will be set to 0 (NA does not work with barplot)
      superabundances=matrix(0,nrow=length(unique(groups))*(topTaxa+1),ncol=ncol(abundances))
      if(aggregate != "none"){
        superabundances=matrix(0,nrow=length(unique(groups))*(topTaxa+1),ncol=length(unique(groups)))
      }
      superrownames=c()
      # loop groups
      prevRows=0
      prevCols=0
      if(aggregate != "none"){
        prevCols=1
      }
      for(group in unique(groups)){
        group.member.indices=which(groups==group)
        #print(paste("Group-wise summing: Number of samples in group",group,": ",length(group.member.indices)))
        sorted.group=sortTaxa(abundances[,group.member.indices], topTaxa = topTaxa)
        aggregated=c()
        superrownames=c(superrownames,rownames(sorted.group))
        if(aggregate != "none"){
          if(aggregate=="mean"){
            aggregated=apply(sorted.group,1,mean) # compute mean row-wise
          }else if(aggregate=="median"){
            aggregated=apply(sorted.group,1,median) # compute median row-wise
          }
          print(paste("Top taxa in group", group))
          print(aggregated)
          superabundances[(prevRows+1):(prevRows+topTaxa+1),prevCols]=aggregated
          prevCols=prevCols+1
        }else{
          superabundances[(prevRows+1):(prevRows+topTaxa+1),(prevCols+1):(prevCols+ncol(sorted.group))]=sorted.group
          prevRows=prevRows+topTaxa+1
          prevCols=prevCols+ncol(sorted.group)
        }
      }
      rownames(superabundances)=superrownames
      if(aggregate != "none"){
        colnames(superabundances)=unique(groups)
      }else{
        colnames(superabundances)=colnames(abundances)
      }
  }else{
    sorted=sortTaxa(abundances,topTaxa=topTaxa)
  }

  if(length(hide.taxa)>0){
    other.index=which(rownames(sorted)=="Others")
    if(length(other.index)>0){
      #print(dim(hidden.taxa))
      hidden.sum=colSums(hidden.taxa)
      sorted[other.index,]=sorted[other.index,]+hidden.sum
    }else{
      warning("Could not add hidden taxa to other taxa!")
    }
  }

  if(sumGroupwise){
    sorted=superabundances
    color.res=getColorVectorGivenTaxaAndColorMap(rownames(sorted),taxon.color.map = taxon.color.map, summedTaxonColor = summedTaxonColor, extendColorMap = extendTaxonColorMap, sumGroupwise = TRUE)
    taxon.color.map=color.res$colormap
    taxon.colors=color.res$colors
  }else{
    color.res=getColorVectorGivenTaxaAndColorMap(rownames(sorted),taxon.color.map = taxon.color.map, summedTaxonColor = summedTaxonColor, extendColorMap = extendTaxonColorMap)
    taxon.color.map=color.res$colormap
    taxon.colors=color.res$colors
  }


  # add margin for the legend
  maxchars=0
  for(taxon.name in rownames(sorted)){
    chars=nchar(taxon.name)
    if(chars>maxchars){
      maxchars=chars
    }
  }
  updated.mar=prev.mar
  #print(paste("Maximal character number",maxchars))
  updated.mar[4]=maxchars*mar.scale
  par(mar=updated.mar,srt=90, las=2)
  #print(updated.mar)

  if(groupNum>1){
    group.colors=assignColorsToGroups(groups,myColors = group.color.map)
    #print(unique(group.colors))
  }else{
    group.colors="black"
  }

  # do the bar plot
  # note that ylim removes the margin definition, so is omitted
  # check presence of ylab and add a default if absent
  if(!("ylab" %in% names(match.call(expand.dots=TRUE)))){
    midpoints=barplot(sorted,col=taxon.colors,xaxt='n',ylab="Abundance",cex.names=0.8, cex.axis=0.8, ...)
  }else{
    midpoints=barplot(sorted,col=taxon.colors,xaxt='n',cex.names=0.8, cex.axis=0.8, ...)
  }
  #print(colnames(sorted))
  mtext(colnames(sorted),col=group.colors,side=1, cex=0.8, line=0.5, at=midpoints)
  prev.xpd=par()$xpd
  par(las=1,srt=0,mar=updated.mar, xpd=NA)

  # add the legend
  if(legend==TRUE){
    # found in: https://stackoverflow.com/questions/42075751/calculating-the-appropriate-inset-value-for-legends-automatically
    coord=par("usr")
    if(sumGroupwise){
      legend(x=coord[2]*legend.shift,y=coord[4],legend=unique(rownames(sorted)),cex=0.8, bg = "white", text.col=unique(taxon.colors))
    }else{
      legend(x=coord[2]*legend.shift,y=coord[4],legend=rownames(sorted),cex=0.8, bg = "white", text.col=taxon.colors)
    }
    #legend(x="topleft",inset=c(legend.shift,0),legend=rownames(sorted),cex=0.8, bg = "white", text.col=taxon.colors)
    #legend(updated.mar[4]*(legend.shift/mar.scale),1,legend=rownames(sorted),cex=0.8, bg = "white", text.col=taxon.colors,xpd=TRUE)
    if(groupNum>1 && legend.hidegroups==FALSE){
      legend(x="bottomleft",inset=c(legend.shift,0),legend=unique(groups),cex=0.8,bg="white", text.col=unique(group.colors))
      #legend(updated.mar[4]*(legend.shift/mar.scale),0,legend=unique(groups),cex=0.8,bg="white", text.col=unique(group.colors),xpd=TRUE)
    }
  }
  # restore par to default values
  par(mar=prev.mar, xpd=prev.xpd)

  if(extendTaxonColorMap){
    return(taxon.color.map)
  }
}

# sort taxa in abundance table xgroup by summed abundance and sum all non-topTaxa into a single group
sortTaxa<-function(xgroup, topTaxa=0){
  # sort taxa by abundance within group
  rowsums=apply(xgroup,1,sum)
  sorted=sort(rowsums,decreasing=TRUE,index.return=TRUE)
  sub.xgroup=xgroup[sorted$ix[1:topTaxa],]
  # sum remaining taxa
  misc=xgroup[sorted$ix[(topTaxa+1):nrow(xgroup)],]
  if(topTaxa<nrow(xgroup)){
    if(ncol(xgroup) > 1){
      misc.summed=apply(misc,2,sum)
      sub.xgroup=rbind(sub.xgroup,misc.summed)
    }else{
      misc.summed=sum(misc)
      sub.xgroup=as.matrix(c(sub.xgroup,misc.summed))
      rownames(sub.xgroup)=c(rownames(xgroup)[sorted$ix[1:topTaxa]],"")
    }
  }
  if(topTaxa<nrow(xgroup)){
    rownames(sub.xgroup)[nrow(sub.xgroup)]="Others"
  }
  sub.xgroup=as.matrix(sub.xgroup)
  return(sub.xgroup)
}

# return a color object with a vector of colors as first entry and a color map as second entry
# given a list of taxon names and optionally a predefined color map
getColorVectorGivenTaxaAndColorMap<-function(namesTopTaxa=c(),taxon.color.map=list(), summedTaxonColor="#a9a9a9", extendColorMap=FALSE, sumGroupwise=FALSE){
  res=list()
  colorMapProvided=TRUE
  # check if taxon color map was provided and initialise one if not
  if(is.null(taxon.color.map) || length(taxon.color.map)==0){
    colorMapProvided=FALSE
    taxon.color.map=list()
  }
  colornumber=length(namesTopTaxa)
  if(sumGroupwise){
    # determine how many non-overlapping top taxa there are
    colornumber=length(unique(namesTopTaxa))
    if("Others" %in% namesTopTaxa){
      colornumber=colornumber-1
    }
  }
  col.vec = seq(0,1,1/colornumber)
  my.colors = hsv(col.vec)
  colors=c()
  color.index=1

  if(sumGroupwise && !colorMapProvided){
    # a color map is needed to deal with multiple assignments
    color.counter=1
    for(name in namesTopTaxa){
      if(!(name %in% names(taxon.color.map)) && !(name=="Others")){
        taxon.color.map[[name]]=my.colors[color.counter]
        color.counter=color.counter+1
      }
    }
    colorMapProvided=TRUE
  } # deal with overlapping names

  if(extendColorMap){
    # generate additional colors if needed
    if(colorMapProvided){
      # number of top taxa with a color assigned
      taxaUnassigned=setdiff(namesTopTaxa,names(taxon.color.map))
      # remove all colors already in custom color map
      colorsAssigned=c(unique(unlist(taxon.color.map)),summedTaxonColor)
      colorsUnassigned=setdiff(unique(my.colors),colorsAssigned)
      increment=5
      # come up with more colors until we have enough
      while(length(colorsUnassigned)<length(taxaUnassigned)){
        col.vec = seq(0,1,1/(colornumber+increment))
        my.colors = hsv(col.vec)
        colorsUnassigned=setdiff(unique(my.colors),colorsAssigned)
        increment=increment+5
      }
      my.colors=colorsUnassigned
    } # end add colors
  }

  # if not yet there, add color for summed taxa
  if(!(summedTaxonColor %in% unlist(taxon.color.map))){
    # gray color for non-top taxa
    taxon.color.map[["Others"]]=summedTaxonColor
  }

  #print(taxon.color.map)

  # loop top taxa
  for(name in namesTopTaxa){
    #print(name)
    # no taxon color map provided: color is taken from color vector
    if(!colorMapProvided){
      if(name!="Others"){
        colors=c(colors, my.colors[color.index])
        taxon.color.map[[name]]=my.colors[color.index]
        color.index=color.index+1
      }else{
        colors=c(colors, summedTaxonColor)
      }
    }
    # color map provided
    else{
      # taxon found in color map
      if(name %in% names(taxon.color.map)){
        #print(paste("Found color",taxon.color.map[[name]],"for name",name))
        colors=c(colors,taxon.color.map[[name]])
      }else{
        # taxon not found, but extend is true
        if(extendColorMap){
          if(name=="Others"){
            taxon.color.map[[name]]=summedTaxonColor
            colors=c(colors, summedTaxonColor)
          }else{
            taxon.color.map[[name]]=my.colors[color.index]
            colors=c(colors, my.colors[color.index])
            color.index=color.index+1
          }
        }else{
          # taxon absent from color map
          print(paste("Taxon",name,"is not present in the color map. It will be assigned the color",summedTaxonColor))
          taxon.color.map[[name]]=summedTaxonColor
          colors=c(colors, summedTaxonColor)
        }
      } # taxon absent from color map
    } # color map provided
  }# loop taxa
  #print(namesTopTaxa)
  #print(colors)
  res[["colors"]]=colors
  res[["colormap"]]=taxon.color.map
  return(res)
}
