#' @title Compare the distribution of a taxon across two groups
#'
#' @description Plot taxon abundances as histograms and calculate the
#' unpaired two-sided Wilcoxon p-value.
#'
#' @param x a matrix with taxa as rows and samples as columns
#' @param taxon the index or name of a taxon to be compared
#' @param groups a vector that provides for each sample its group membership (numeric or character)
#' @param group1 the identifier of the first group (numeric or character)
#' @param group2 the identifier of the second group (numeric or character)
#' @param name1 the name to be displayed for the first group
#' @param name2 the name to be displayed for the second group
#' @examples
#' data("ibd_taxa")
#' data("ibd_metadata")
#' groups=as.vector(ibd_metadata$Diagnosis)
#' taxon="Faecalibacterium_prausnitzii"
#' compareDistribs(ibd_taxa,taxon=taxon,groups=groups,group1="UC",group2="Control")
#' groups[groups=="UC"]="IBD"
#' groups[groups=="CD"]="IBD"
#' # Faecalibacterium has a significantly higher relative abundance in the control group
#' compareDistribs(ibd_taxa,taxon=taxon,groups=groups,group1="IBD",group2="Control")
#' @export
compareDistribs<-function(x, taxon, groups=c(), group1=1, group2=2, name1=as.character(group1), name2=as.character(group2)){
  indices.group1=which(groups==group1)
  indices.group2=which(groups==group2)
  index=NA
  if(!is.numeric(taxon)){
    index=which(rownames(x)==taxon)
  }else{
    index=taxon
  }
  #print(index)
  abundances1=as.numeric(x[index,indices.group1])
  abundances2=as.numeric(x[index,indices.group2])
  w.out=wilcox.test(abundances1,abundances2)
  title=paste(rownames(x)[index],", \nWilcoxon p-value:",round(w.out$p.value,3),sep="")
  compareDistribsPure(values1=abundances1, values2=abundances2, name1=name1,name2=name2, main=title, displayW=FALSE)
}

# values1: first distribution
# values2: second distribution
# name1: name of first distribution
# name2: name of second distribution
# xlab: x axis label
# legend.position: position of the legend
# breaks: breaks to use for hist
# displayW: append Wilcoxon p-value to main
compareDistribsPure<-function(values1,values2, main="Comparison", name1="", name2="", xlab="Abundance", legend.position="topright", breaks="FD", displayW=TRUE){
  w.out=wilcox.test(values1,values2)
  col1=rgb(0,1,0,0.5)
  col2=rgb(1,0,0,0.5)
  # limits
  xmax=max(values1,na.rm=TRUE)
  xmin=min(values1,na.rm=TRUE)
  ymax=max(values2,na.rm=TRUE)
  ymin=min(values2,na.rm=TRUE)
  max=max(xmax,ymax)
  min=min(ymin,xmin)
  xmaxD=max(hist(values1,breaks=breaks,plot=FALSE)$density)
  ymaxD=max(hist(values2,breaks=breaks,plot=FALSE)$density)
  maxD=max(xmaxD,ymaxD)
  #maxD=maxD+0.5 # add a margin
  if(displayW){
    if(nzchar(main)){
      title=paste(main,", \nWilcoxon p-value:",round(w.out$p.value,3),sep="")
    }else{
      title=paste("Wilcoxon p-value:",round(w.out$p.value,3),sep="")
    }
  }else{
    title=main
  }
  hist(values1,breaks=breaks,xlim=c(min,max), ylim=c(0,maxD), prob=TRUE,col=col1, border=col1,xlab=xlab, main=title)
  hist(values2,breaks=breaks,prob=TRUE,col=col2, border=col2,add=TRUE)
  legend(legend.position,legend=c(name1,name2), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")
}
