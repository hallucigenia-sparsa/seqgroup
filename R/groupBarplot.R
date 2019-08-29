#' @title Bar plot of taxon composition with group support
#' @description Sort taxa by summed abundance across samples and plot sorted taxon composition with a bar per sample
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param groups group membership vector with as many entries as samples
#' @param taxon.color.map map of taxon-specific colors, should match row names; taxa not present in the color map will be colored in summedTaxonColor
#' @param group.color.map map of group-specific colors, should match group names
#' @param topTaxa number of top taxa to be plotted
#' @param sortGroupwise if true, samples are sorted according to groups
#' @param randSampleNum if larger 0, sortGroupwise is set to true and the indicated sample number is randomly selected for each group
#' @param summedTaxonColor the color of the summed taxa, by default gray
#' @param extendTaxonColorMap if true, taxa not in the taxon color map are added there and the extended color map is returned
#' @param legend add a legend with the color code
#' @param legend.shift increase/decrease this parameter to shift the color legend further to the right/left
#' @param \\dots Additional arguments passed to plot()
#' @param return if extendTaxonColorMap is true, the taxon color map is returned
#' @examples
#' data(ibd_taxa)
#' data(ibd_metadata)
#' # taxon abundances were prefiltered and therefore do not add up to 1
#' groupBarplot(ibd_taxa,groups=as.vector(ibd_metadata$Diagnosis),randSampleNum=7)
#' data(ibd_lineages)
#' ibd_genera=aggregateTaxa(ibd_taxa,ibd_lineages,taxon.level = "genus")
#' groupBarplot(ibd_genera,groups=as.vector(ibd_metadata$Diagnosis),randSampleNum=7)
#' @export
groupBarplot<- function(abundances, groups=c(), taxon.color.map=NULL, group.color.map=NULL, topTaxa=10, sortGroupwise=TRUE, randSampleNum=NA, summedTaxonColor="#a9a9a9", extendTaxonColorMap=FALSE, legend=TRUE, legend.shift=1, ...){

  if(length(groups)>0){
    if(length(groups)!=ncol(abundances)){
      stop("Each sample should have a group assigned.")
    }
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

  groupNum=length(unique(groups))
  prev.mar=par()$mar
  mar.scale=0.5 # maximal number of characters is scaled by this number to compute mar on the right side

  sorted=NULL
  taxon.colors=c()
  group.colors=c()

  # sort samples according to group membership
  # TODO: only collecting indices is sufficient, no need to use cbind
  if(sortGroupwise && (groupNum>1 || randSampleNum>0)){
    counter=1
    abundancesSorted=NULL
    updated.groups=c()
    # loop groups
    for(group in unique(groups)){
      group.member.indices=which(groups==group)
      print(paste("Number of samples in group",group,":",length(group.member.indices)))
      if(randSampleNum>0){
        if(length(group.member.indices)>randSampleNum){
          group.member.indices=sample(group.member.indices)[1:randSampleNum]
        }
      }
      xgroup=as.matrix(abundances[,group.member.indices])
      if(length(group.member.indices)==1){
        colnames(xgroup)=colnames(abundances)[group.member.indices]
      }
      #print(colnames(xgroup))
      updated.groups=c(updated.groups,rep(group,length(group.member.indices)))
      if(counter>1){
        abundancesSorted=cbind(abundancesSorted,xgroup)
      }else{
        abundancesSorted=xgroup
      }
      counter=1+counter
    } # end groups loop
    abundances=abundancesSorted
    groups=updated.groups
    #print(dim(abundances))
    #print(groups)
  }

  sorted=sortTaxa(abundances,topTaxa=topTaxa)
  color.res=getColorVectorGivenTaxaAndColorMap(rownames(sorted),taxon.color.map = taxon.color.map, summedTaxonColor = summedTaxonColor, extendColorMap = extendTaxonColorMap)
  taxon.color.map=color.res$colormap
  taxon.colors=color.res$colors

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
  # par(las=2, srt=90, mar = c(5, 5, 4, 4))

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
    midpoints=barplot(sorted,col=taxon.colors,xaxt='n',ylab="Abundance",cex.names=0.8, ...)
  }else{
    midpoints=barplot(sorted,col=taxon.colors,xaxt='n',cex.names=0.8, ...)
  }
  #print(colnames(sorted))
  mtext(colnames(sorted),col=group.colors,side=1, cex=0.8, line=0.5, at=midpoints)
  par(las=1,srt=0)

  # add the legend
  if(legend==TRUE){
    legend(updated.mar[4]*(legend.shift/mar.scale),1,legend=rownames(sorted),cex=0.8, bg = "white", text.col=taxon.colors,xpd=TRUE)
    if(groupNum>1){
      legend(updated.mar[4]*(legend.shift/mar.scale),0,legend=unique(groups),cex=0.8,bg="white", text.col=unique(group.colors),xpd=TRUE)
    }
  }
  par(mar=prev.mar)

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
getColorVectorGivenTaxaAndColorMap<-function(namesTopTaxa=c(),taxon.color.map=list(), summedTaxonColor="#a9a9a9", extendColorMap=FALSE){
  res=list()
  colorMapProvided=TRUE
  # check if taxon color map was provided and initialise one if not
  if(is.null(taxon.color.map) || length(taxon.color.map)==0){
    colorMapProvided=FALSE
    taxon.color.map=list()
  }
  colornumber=length(namesTopTaxa)
  col.vec = seq(0,1,1/colornumber)
  my.colors = hsv(col.vec)
  colors=c()
  color.index=1

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
  # loop top taxa
  for(name in namesTopTaxa){
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
        #print(paste("Found color",colormap[[name]],"for name",name))
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
  res[["colors"]]=colors
  res[["colormap"]]=taxon.color.map
  return(res)
}
