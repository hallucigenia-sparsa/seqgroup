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
#'
compareDistribs<-function(x, taxon, groups=c(), group1=1, group2=2, name1="group1", name2="group2"){
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
  col1=rgb(0,1,0,0.5)
  col2=rgb(1,0,0,0.5)
  # limits
  xmax=max(abundances1,na.rm=TRUE)
  xmin=min(abundances1,na.rm=TRUE)
  ymax=max(abundances2,na.rm=TRUE)
  ymin=min(abundances2,na.rm=TRUE)
  max=max(xmax,ymax)
  min=min(ymin,xmin)
  xmaxD=max(hist(abundances1,breaks="FD",plot=FALSE)$density)
  ymaxD=max(hist(abundances2,breaks="FD",plot=FALSE)$density)
  maxD=max(xmaxD,ymaxD)
  #maxD=maxD+0.5 # add a margin
  title=paste(rownames(x)[index],", \nWilcoxon p-value:",round(w.out$p.value,3),sep="")
  hist(abundances1,breaks="FD",xlim=c(min,max), ylim=c(0,maxD), prob=TRUE,col=col1, border=col1,xlab="Abundance", main=title)
  hist(abundances2,breaks="FD",prob=TRUE,col=col2, border=col2,add=TRUE)
  legend("right",legend=c(name1,name2), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")
}