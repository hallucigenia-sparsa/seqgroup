#' @title PCoA for microbial sequencing data
#'
#' @description A wrapper around various PCoA-based analyses implemented in vegan. The wrapper can handle groups and
#' metadata. The na.action is set to na.omit, however envfit cannot deal with
#' missing values, therefore if metadata are provided, they should be free of missing values.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param metadata an optional data frame with metadata items, where samples are in the same order as in x, if provided and rda is FALSE, envfit is carried out
#' @param groupAttrib optional: the name of a metadata item that refers to a vector that provides for each sample its group membership
#' @param groups an optional vector that provides for each sample its group membership and which is overridden by groupAttrib, if provided
#' @param clusters an optional vector that provides for each sample its cluster membership as an integer (cluster membership is visualized through shape, up to 10 different shapes are possible)
#' @param labels an optional vector that provides for each sample a label to display
#' @param time an optional vector with as many time points as samples, adds arrows between consecutive time points (time points should be in ascending order)
#' @param hiddenSamples an optional vector with indices of samples to be hidden (they will be taken into account for PCoA/RDA/envfit, but are not displayed)
#' @param dis dissimilarity or distance supported by vegan's vegdist function
#' @param rda carry out an RDA instead of a PCoA
#' @param scale scale numeric metadata (subtract the mean and divide by standard deviation)
#' @param doScree do a Scree plot
#' @param topTaxa if larger than zero: show the top N taxa most correlated to principal components as arrows in the PCoA
#' @param topMetadata if larger than zero, metadata provided and rda false: show the top N most significant numeric metadata as arrows and the top N most significant factor metadata as text in the PCoA
#' @param arrowFactor the length of taxon arrows (determined by scaled covariance) is multiplied with this factor
#' @param metadataFactor the length of numeric metadata arrows (determined by Pearson correlation) is multiplied with this factor
#' @param centroidFactor centroid positions are multiplied with this factor
#' @param taxonColor the color of the taxon arrows and text
#' @param metadataColor the color of the metadata arrows and text
#' @param xlim range shown on the x axis
#' @param ylim range shown on the y axis
#' @param permut number of permutations in envfit
#' @param pAdjMethod method for multiple testing correction supported by p.adjust for envfit p-values
#' @param qvalThreshold threshold on multiple-testing corrected envfit p-values
#' @param dimensions the principal components used for plotting, by default the first and second
#' @param \\dots Additional arguments passed to plot()
#' @export
#'

seqPCoA<-function(abundances, metadata=NULL, groupAttrib="", groups=c(), clusters=c(), labels=c(), time=c(), hiddenSamples=c(), dis="bray", rda=FALSE, scale=FALSE, doScree=FALSE, topTaxa=10, topMetadata=10, arrowFactor=0.5, metadataFactor=1, centroidFactor=1, taxonColor="brown", metadataColor="blue", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3), permut=1000, pAdjMethod="BH", qvalThreshold=0.05, dimensions=c(1,2), ...){

  if(rda && is.null(metadata)){
    stop("Metadata are needed for RDA!")
  }

  if(!is.null(metadata) && scale){
    numeric.metadata=getMetadataSubset(metadata,type="numeric")
    catbin.metadata=getMetadataSubset(metadata,type="catbin")
    scaled.numeric.metadata=scale(numeric.metadata,center=TRUE, scale=TRUE)
    metadata=cbind(scaled.numeric.metadata,catbin.metadata)
  }

  pch.value=16

  if(rda){
    pcoa.res=capscale(data.frame(t(abundances))~.,metadata,distance=dis, na.action = "na.omit")

  }else{

    # carry out PCoA
    pcoa.res=capscale(data.frame(t(abundances))~1,distance=dis, na.action = "na.omit")
    if(!is.null(metadata) && topMetadata>0){

      # carry out envfit
      ef=envfit(pcoa.res,metadata,perm=permut)

      # correct for multiple testing using code in http://www.davidzeleny.net/anadat-r/doku.php/en:indirect_ordination_suppl
      pvals.vectors=p.adjust(ef$vectors$pvals, method=pAdjMethod)
      pvals.factors=p.adjust(ef$factors$pvals, method=pAdjMethod)
      # sort p-values in ascending orders, keep sorting result for other vector output
      pvals.vectors.sorted=sort(pvals.vectors, index.return=TRUE)
      pvals.vectors=pvals.vectors.sorted$x
      # factors make only use of names, but not of any other result
      pvals.factors=sort(pvals.factors)
      indices.vectors=which(pvals.vectors<qvalThreshold)
      indices.factors=which(pvals.factors<qvalThreshold)
      print(paste(length(indices.vectors),"significant numeric metadata found, in order of significance:"))
      #print(pvals.vectors)
      for(index in indices.vectors){
        #print(names(ef$vectors$r[index]))
        print(names(ef$vectors$r[pvals.vectors.sorted$ix[index]]))
      }
      print(paste(length(indices.factors),"significant categoric metadata found, in order of significance:"))
      for(index in indices.factors){
        print(names(pvals.factors)[index])
      }
      if(length(indices.vectors)>topMetadata){
        indices.vectors=indices.vectors[1:topMetadata]
      }
      if(length(indices.factors)>topMetadata){
        indices.factors=indices.factors[1:topMetadata]
      }
    } # do envfit
  } # not rda

  if(doScree){
    barplot(pcoa.res$CA$eig, names.arg=c(1:length(pcoa.res$CA$eig)), main="Scree plot", xlab="Eigenvalue index", ylab="Eigenvalue")
  }

  # assign the right labels
  xlab=paste("PCoA",dimensions[1]," [",round(pcoa.res$CA$eig[dimensions[1]],2),"]",sep="")
  ylab=paste("PCoA",dimensions[2]," [",round(pcoa.res$CA$eig[dimensions[2]],2),"]",sep="")

  colors="red"
  if(groupAttrib!=""){
    groups=metadata[[groupAttrib]]
    colors=assignColorsToGroups(groups)
  }else if(length(groups)>0){
    colors=assignColorsToGroups(groups)
  }

  # select indices of samples to show (if hiddenSamples is empty, this will be all)
  selected.sample.indices=setdiff(1:ncol(abundances),hiddenSamples)

  # 15=square, 16=circle, 17=triangle point up, 18=diamond, 25=triangle point down,
  # 3=plus sign, 4=multiplier sign, 8=star sign, 9=diamond with plus sign, 7=square with plus sign
  clus.pch.values=c(15,16,17,18,25,3,4,8,9,7)
  # assign cluster memberships as shapes
  if(length(clusters)>0){
    clusnum=length(unique(clusters))
    if(clusnum>length(clus.pch.values)){
      stop(paste("No more than",length(clus.pch.values),"clusters can be visualized."))
    }
    pch.value=c()
    for(clus.index in 1:length(clusters)){
      pch.prop=clus.pch.values[clusters[clus.index]]
      pch.value=c(pch.value,as.numeric(pch.prop))
    }
    pch.value=pch.value[selected.sample.indices]
    #print(pch.value[20:50])
  }

  plot(pcoa.res$CA$u[selected.sample.indices,dimensions], col=colors[selected.sample.indices], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, bg=colors[selected.sample.indices], pch=pch.value, ...)

  # add sample labels if requested
  if(length(labels)>0){
    text(pcoa.res$CA$u[selected.sample.indices,dimensions[1]],pcoa.res$CA$u[selected.sample.indices,dimensions[2]],labels=labels[selected.sample.indices], pos=3, cex=0.9)
  }

  # add arrows between consecutive time points if requested
  if(length(time)>0){
    for(i in 1:(nrow(pcoa.res$CA$u)-1)){
      # both samples to be connected are visible
      if(!(i %in% hiddenSamples) && !((i+1) %in% hiddenSamples)){
        if(length(groups)==0 || groups[i]==groups[i+1]){
          arrows(x0=pcoa.res$CA$u[i,dimensions[1]],y0=pcoa.res$CA$u[i,dimensions[2]],x1=pcoa.res$CA$u[i+1,dimensions[1]],y1=pcoa.res$CA$u[i+1,dimensions[2]], length=0.1, col="gray", lty=2, angle=20) # length=length of the edges of the arrow head in inches, lty=2 is dashed, angle refers to shaft vs edge of arrow head
        } # within the same group
      } # visible
    } # loop samples
  } # time provided

  if(topTaxa>0){
    # taken from biplot.pcoa in ape
    # columns are taxa
    Y=t(abundances)
    k <- ncol(pcoa.res$CA$u)
    n <- nrow(Y)
    # standardize eigen vectors
    points.stand <- scale(pcoa.res$CA$u[,dimensions])
    # covariance between taxa and selected principal components (standardized eigen vectors)
    S <- cov(Y, points.stand)
    # scale S by the eigen values
    U <- S %*% diag((pcoa.res$CA$eig[dimensions]/(n-1))^(-0.5))
    colnames(U) <- colnames(pcoa.res$CA$u[,dimensions])

    # select top taxa in U
    # arrrow length in either direction codes strength of correlation to eigen vectors
    rowsums=apply(abs(U),1,sum)
    sorted=sort(rowsums,index.return=TRUE,decreasing=TRUE)
    U.selected=U[sorted$ix[1:topTaxa],]

    arrows(0, 0, U.selected[, 1] * arrowFactor, U.selected[, 2] * arrowFactor, col = taxonColor,length = 0.1, lty=2)

    shift=0.1
    for(i in 1:nrow(U.selected)){
      for(j in 1:ncol(U.selected)){
        if(U.selected[i,j]>0){
          U.selected[i,j]=U.selected[i,j]+shift
        }else if(U.selected[i,j]<0){
          U.selected[i,j]=U.selected[i,j]-shift
        }
      }
    }
    text(U.selected*arrowFactor, rownames(U.selected), cex = 0.9, col = taxonColor)
  }

  if(!is.null(metadata) && topMetadata>0){
    if(length(indices.vectors)>0){
      # add arrows for numeric metadata
      ef.arrows=as.matrix(ef$vectors$arrows[pvals.vectors.sorted$ix[indices.vectors],dimensions])
      if(length(indices.vectors)==1){
        ef.arrows=t(ef.arrows)
      }
      lengths=ef$vectors$r[pvals.vectors.sorted$ix[indices.vectors]]
      ef.names=names(ef$vectors$r[pvals.vectors.sorted$ix[indices.vectors]])
      arrows(0, 0, ef.arrows[, 1]*lengths*metadataFactor, ef.arrows[, 2]*lengths*metadataFactor, col = metadataColor, length=0.1)

      # move metadata labels to nicer places
      shift=0.05
      for(i in 1:nrow(ef.arrows)){
        for(j in 1:ncol(ef.arrows)){
          if(ef.arrows[i,j]>0){
            ef.arrows[i,j]=ef.arrows[i,j]+shift
          }else if(ef.arrows[i,j]<0){
            ef.arrows[i,j]=ef.arrows[i,j]-shift
          }
        }
      }
      text(ef.arrows*lengths*metadataFactor, ef.names, cex = 0.8, col = metadataColor)
    }
    if(length(indices.factors)>0){
      # the order of factors was altered, but since we match by name, this is not problematic
      sig.factors=names(indices.factors)
      # find centroids belonging to significant factors
      # centroid: mean or median dissimilarity of samples belonging to given level of a factor in PCoA
      # significance of centroid separation: TukeyHSD?
      for(sig.factor in sig.factors){
        sig.factor.centroid.indices=which(ef$factors$var.id==sig.factor)
        # several values of the categoric variable can be significant
        for(sig.factor.centroid.index in sig.factor.centroid.indices){
          sig.factor.value=rownames(ef$factors$centroids)[sig.factor.centroid.index]
          centroid.pos=ef$factors$centroids[sig.factor.centroid.index,]
          text(x=centroid.pos[1]*centroidFactor,y=centroid.pos[2]*centroidFactor,sig.factor.value, cex=0.8,col=metadataColor)
        } # end loop values of significant factors
      } # end loop significant factors
    } # significant factors found
  }

  if(length(groups)>0){
    legend("topright",legend=unique(groups[selected.sample.indices]),cex=0.9, pch = rep("*",length(unique(groups[selected.sample.indices]))), col = unique(colors[selected.sample.indices]), bg = "white", text.col="black")
  }
  if(length(clusters)>0){
    legend("topleft",legend=unique(clusters[selected.sample.indices]),cex=0.9, pch = unique(pch.value), col = "black", bg = "white", text.col="black")
  }

}

# expects group membership vector as input and returns a color vector
# assign the same color to members of the same group
# color vector is as long as group membership vector
# if returnMap is true, the color map is returned
assignColorsToGroups<-function(groups, my.color.map = list(), returnMap=FALSE){
  groupNum=length(unique(groups))
  col.vec = seq(0,1,1/groupNum)
  hues = hsv(col.vec)
  #print(hues)
  hueCounter=1
  colors=c()
  # fill the color map
  if(length(my.color.map)==0){
    for(group.index in 1:length(groups)){
      group=as.character(groups[group.index])
      if(!(group %in% names(my.color.map))){
        my.color.map[[group]]=hues[hueCounter]
        hueCounter=hueCounter+1
      }
    } # loop samples
  } # color map already filled

  # assign colors from the color map
  for(group.index in 1:length(groups)){
    colors=c(colors,my.color.map[[groups[group.index]]])
  }
  if(returnMap){
    res=list(colors,my.color.map)
    names(res)=c("colors","colorMap")
    return(res)
  }
  return(colors)
}


