#' @title PCoA for microbial sequencing data
#'
#' @description A wrapper around various PCoA-based analyses implemented in vegan. The wrapper can handle groups and
#' metadata. PCoA is carried out sample-wise. The na.action is set to na.omit, however envfit cannot deal with
#' missing values, therefore if metadata are provided, they should be free of missing values.
#'
#' @details When a reference and groups are provided and the number of group memberships does not equal the number of samples in the combined abundance table,
#' groups are automatically extended such that reference samples are assigned to a single group with name refName, which is colored in gray.
#' The color vector is likewise extended if provided. If a clusters, time and/or labels vector is provided together with a referene, it has to refer to both data sets.
#' Samples in the abundance matrix are appended after the reference samples, so cluster memberships, time points and/or labels have to be provided in the same order.
#' Different total counts in abundances and reference samples may bias the result, therefore rarefyRef allows to rarefy both to the same total count after matching.
#' RarefyRef will rarefy all samples to the lowest total count found in any sample. Rows with zero counts after rarefaction are removed.
#' When topTaxa is set larger zero, significant top-varying taxa are shown. The permutation test is carried out by shuffling the selected number of top-covarying taxa.
#' Multiple testing correction on parameter-free p-values is then only applied to these top-covarying taxa. The strength of covariance is determined as
#' the norm of the vectors resulting from multiplying the standardized eigen vectors with the taxa.
#' In contrast, envfit p-values are computed for all metadata and multiple-testing correction is consequently applied to all metadata provided, though
#' only the selected number of most significant metadata are shown. Thus, topTaxa ranks taxa by covariance with significance as a filter, whereas
#' topMetadata ranks metadata by significance.
#' The number in the axis label brackets refers to the proportion of variance explained as computed with vegan's eigenvals function.
#' Note that ordination values in plots are not scaled (equivalent to plot.cca(scaling=0)).
#' If groups are provided, drawEllipse can be enabled to draw ellipses with vegan's ordiellipse and to run vegan's adonis (PERMANOVA) to assess whether group composition
#' differs significantly. Adonis is sensitive to differences in within-group variation. To assess whether within-group variation differs, betadisper can be carried out by
#' enabling groupDispersion.
#' In addition, if groups are provided, seqPCoA computes cluster quality indices to assess group separation with package clusterCrit.
#' Concerning interpretation of cluster quality indices: The silhouette index ranges between -1 and 1, a large Dunn's, Silhouette or Calinski-Harabasz index and a small
#' Davies-Bouldin or C-index indicate well-defined clusters, respectively.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param reference an optional reference data set on which abundances are mapped; data are merged by matching row names (nonmatching ones are kept as sum); cannot be combined with rda or topMetadata (topMetadata needs to be set to zero)
#' @param rarefyRef rarefy abundance and reference samples to the minimum total count found in any of the samples; recommended when the total counts differ
#' @param addToRefStepwise compute ordination coordinates for each sample added to the reference separately; cannot be combined with metadata, rda or drawEllipse
#' @param refName group name for reference samples
#' @param metadata an optional data frame with metadata items as columns, where samples are in the same order as in abundances and data types (factor vs numeric) are supposed to be correct; if provided and rda is FALSE, envfit is carried out
#' @param groupAttrib optional: the name of a metadata item that refers to a vector that provides for each sample its group membership
#' @param groups an optional vector that provides for each sample its group membership and which is overridden by groupAttrib, if provided
#' @param groupColors an optional map of predefined colors for groups that matches names in groups (which should be strings); if reference is provided, refName is added if absent
#' @param colors an optional vector of colors to be used to color samples; it overrides groupColors if provided
#' @param clusters an optional vector that provides for each sample its cluster membership (cluster membership is visualized through shape, up to 10 different shapes are possible)
#' @param labels an optional vector that provides for each sample a label to display
#' @param sizes a vector of  numeric values that will be displayed as varying sample sizes (sizes will be shifted into positive range if necessary and scaled between 0.5 and 2.5)
#' @param size.legend a string displayed as a legend for size
#' @param time an optional vector with as many time points as samples, adds arrows between consecutive time points (time points should be in ascending order)
#' @param hiddenTaxa an optional vector with names of taxa to be hidden (they will be taken into account for PCoA/RDA/envfit, but are not displayed among top co-varying taxa)
#' @param hiddenSamples an optional vector with indices of samples to be hidden (they will be taken into account for PCoA/RDA/envfit, but are not displayed)
#' @param dis dissimilarity or distance supported by vegan's vegdist function (if set to cor, a PCA is carried out using vegan's function rda with scale set to true)
#' @param rda carry out an RDA instead of a PCoA using vegan's capscale function
#' @param scale scale numeric metadata (subtract the mean and divide by standard deviation)
#' @param doScree do a Scree plot
#' @param topTaxa if larger than zero: show the top N taxa most strongly covarying with principal components as arrows in the PCoA if they are significant according to a permutation test
#' @param topMetadata if larger than zero, metadata provided and rda false: show the top N most significant numeric metadata as arrows and the top N most significant factor metadata as text in the PCoA
#' @param arrowFactor the length of taxon arrows (determined by scaled covariance) is multiplied with this factor
#' @param metadataFactor the length of numeric metadata arrows (determined by Pearson correlation) is multiplied with this factor
#' @param centroidFactor centroid positions (representing categoric metadata) are multiplied with this factor
#' @param taxonColor the color of the taxon arrows and text
#' @param metadataColor the color of the metadata arrows and text
#' @param drawEllipse if groups or groupAttrib given, draw polygons encapsulating groups using vegan's ordiellipse function (kind is sd, conf given via ellipseConf); print adonis R2 and p-value (permutation number given via env.permut)
#' @param clusterQualityIndex if groups or groupAttrib given, report cluster quality according to silhouette function or any criterium supported by package clusterCrit (default: silhouette, set to none to disable computation of cluster quality)
#' @param groupDispersion if groups or groupAttrib given, report Tukey's HSD test on differences in group dispersions (avg distance of group members to centroid) and do boxplot (wraps vegan's betadisper)
#' @param xlim range shown on the x axis, by default the minimum and maximum of the first selected component
#' @param ylim range shown on the y axis, by default the minimum and maximum of the second selected component
#' @param permut number of permutations for top-covarying taxa; if NA, NULL or smaller than 1, no permutation test is carried out
#' @param env.permut number of permutations for envfit, if drawEllipse is true, for adonis
#' @param pAdjMethod method for multiple testing correction supported by p.adjust for top-covarying taxon and envfit p-values
#' @param qvalThreshold threshold on multiple-testing corrected top-covarying taxon and envfit p-values
#' @param ellipseConf confidence limit for drawEllipse
#' @param dimensions the principal components used for plotting, by default the first and second
#' @param \\dots Additional arguments passed to plot()
#' @examples
#' data("ibd_taxa")
#' data("ibd_metadata")
#' ibd_metadata=assignMetadataTypes(ibd_metadata,categoric=c("SRA_metagenome_name","Diagnosis"))
#' seqPCoA(ibd_taxa,groups=as.vector(ibd_metadata$Diagnosis),topTaxa=30, drawEllipse=TRUE)
#' # remove 65 samples with missing calprotectin measurements or other missing values in the metadata
#' na.indices=unique(which(is.na(ibd_metadata),arr.ind=TRUE)[,1])
#' indices.to.keep=setdiff(1:nrow(ibd_metadata),na.indices)
#' ibd_metadata=ibd_metadata[indices.to.keep,]
#' ibd_taxa=ibd_taxa[,indices.to.keep]
#' seqPCoA(ibd_taxa,metadata=ibd_metadata,groups=as.vector(ibd_metadata$Diagnosis),topTaxa=30)
#' @export
#'

seqPCoA<-function(abundances, reference=NULL, rarefyRef=FALSE, addToRefStepwise=FALSE, refName="ref", metadata=NULL, groupAttrib="", groups=c(), groupColors=NULL, colors=c(), clusters=c(), labels=c(), sizes=c(), size.legend="", time=c(), hiddenTaxa=c(), hiddenSamples=c(), dis="bray", rda=FALSE, scale=FALSE, doScree=FALSE, topTaxa=10, topMetadata=10, arrowFactor=0.5, metadataFactor=1, centroidFactor=1, taxonColor="brown", metadataColor="blue", drawEllipse=FALSE, clusterQualityIndex="silhouette", groupDispersion=FALSE, xlim=NULL, ylim=NULL, permut=1000, env.permut=1000, pAdjMethod="BH", qvalThreshold=0.05, ellipseConf=0.95, dimensions=c(1,2), ...){

  # Test
  # path.vdp="/Users/u0097353/Documents/Documents_Karoline/MSysBio_Lab/Results/Nephrology/Data/vdp_genera.txt"
  # vdp=read.table(path.vdp,row.names=1,header=TRUE,sep="\t")
  # data("ibd_lineages")
  # ibd.genera=aggregateTaxa(ibd_taxa, lineages=ibd_lineages,taxon.level="genus")
  # ibd.genera.counts = round(ibd.genera*10000) # scale to counts (vdp already sums to 10K)
  # seqPCoA(ibd.genera.counts,reference=vdp, rarefyRef=TRUE, groups=groups, drawEllipse = TRUE, topMetadata = 0) # 53 matching genera
  # # beautified:
  # seqPCoA(ibd.genera.counts,reference=vdp, rarefyRef=TRUE, groups=groups, drawEllipse = TRUE, topMetadata = 0, xlim=c(-0.1,0.15), ylim=c(-0.1,0.1), arrowFactor=0.0001, clusterQualityIndex = "silhouette")
  # test step-wise
  # seqPCoA(ibd.genera.counts[,1:20],reference=vdp, rarefyRef=TRUE, groups=groups,xlim=c(-0.1,0.15), ylim=c(-0.1,0.1), addToRefStepwise = TRUE)

  if(length(clusters)>0){
    clusters=as.character(clusters)
  }

  metadata.to.plot=c()
  #metadata.to.plot=c("IS_TOTAL","PCSG")

  if(rda && is.null(metadata)){
    stop("Metadata are needed for RDA!")
  }

  if(rda && !is.null(reference)){
    stop("RDA is not supported when a reference is provided.")
  }

  if(addToRefStepwise && is.null(reference)){
    stop("A reference is needed to add samples to reference step-wise.")
  }

  if(doScree && addToRefStepwise){
    stop("Scree plot is not available for step-wise addition of samples to reference.")
  }

  if(rda && addToRefStepwise){
    stop("RDA is not available for step-wise addition of samples to reference.")
  }

  # make it easier for users
  if(addToRefStepwise){
    if(topTaxa>0){
      print("Biplot with taxa is not available for step-wise addition of samples to reference. TopTaxa is set to 0.")
    }
    topTaxa=0
    if(!is.null(metadata)){
      print("Metadata cannot be handled when samples are added step-wise to reference. Metadata are ignored.")
    }
    metadata=NULL
  }

  # if no metadata are provided, we cannot look for top-significant metadata items
  if(is.null(metadata)){
    topMetadata=0
  }

  # TODO: there is no reason why envfit is not supported in case metadata object provides data for both - to be modified
  if(topMetadata > 0 && !is.null(reference)){
    stop("Envfit is not supported when a reference is provided. Please set topMetadata to 0.")
  }

  lastIndexRef=0

  if(!is.null(reference)){
    lastIndexRef=ncol(reference)
    # match abundances and reference by their row names
    res=intersectTables(reference,abundances,byRow = TRUE, keepSumNonMatched = TRUE)
    # append matched abundances to reference
    abundances=cbind(res$table1,res$table2)
    if(rarefyRef){
      # rarefy and discard taxa with zero abundance after rarefaction
      filtered=rarefyFilter(abundances)
      abundances=filtered$rar
      # update lastIndexRef to deal with columns removed during rarefaction
      ref.indices=intersect(1:lastIndexRef, filtered$colindices)
      lastIndexRef=length(ref.indices)
      # check
      print(paste("Name of last reference sample:",colnames(abundances))[lastIndexRef])
    }
    if(length(groups)>0){
      if(length(groups)!=ncol(abundances)){
        # extend groups to reference
        groups=c(rep(refName,ncol(res$table1)),groups)
      }
    }
    if(!is.null(groupColors)){
      if(!(refName %in% names(groupColors))){
        groupColors[[refName]]="gray"
      }
    }
    if(length(colors)>0){
      if(length(colors)!=ncol(abundances)){
        colors=c(rep("gray",ncol(res$table1)),colors)
      }
    }
    if(length(clusters)>0){
      if(length(clusters)!=ncol(abundances)){
        stop("Please provide cluster memberships of reference and abundances combined, with reference first.")
      }
    }
    if(length(labels)>0){
      if(length(labels)!=ncol(abundances)){
        stop("Please provide labels of reference and abundances combined, with reference first.")
      }
    }
    if(length(time)>0){
      if(length(time)!=ncol(abundances)){
        stop("Please provide time points for reference and abundances combined, with reference first.")
      }
    }
  } # reference provided

  if(!is.null(metadata) && scale){
    # print("Scaling metadata")
    numeric.metadata=getMetadataSubset(metadata,type="numeric")
    catbin.metadata=getMetadataSubset(metadata,type="catbin")
    scaled.numeric.metadata=scale(numeric.metadata,center=TRUE, scale=TRUE)
    metadata=cbind(scaled.numeric.metadata,catbin.metadata)
  }

  min.cex=0.5
  max.cex=2.5
  pch.value=16
  # 15=square, 16=circle, 17=triangle point up, 18=diamond, 25=triangle point down,
  # 3=plus sign, 4=multiplier sign, 8=star sign, 9=diamond with plus sign, 7=square with plus sign
  clus.pch.values=c(15,16,17,18,25,3,4,8,9,7)
  display.size.legend=FALSE
  sample.coords=matrix(NA, nrow=ncol(abundances),ncol=2) # samples x axes

  if(rda){

    # it is not possible to give dissimilarities to capscale, so they need to be recomputed in case betadisper is enabled
     pcoa.res=capscale(data.frame(t(abundances))~.,metadata,distance=dis, na.action = "na.omit")

  }else if(addToRefStepwise){

    # loop over samples, carry out PCoA and extract coordinates of joint ordination with reference
    background=abundances[,1:lastIndexRef]
    samples=abundances[,(lastIndexRef+1):ncol(abundances)]
    print(paste("Carrying out",ncol(samples),"ordinations to match samples to reference one by one"))
    # loop samples to match to the background
    for(ord.index in 0:ncol(samples)){
      if(ord.index==0){
        matched=background
      }else{
        matched=cbind(background,samples[,ord.index])
      }
      if(dis=="cor"){
        # carry out standard PCA
        pcoa.res=rda(data.frame(t(matched)), scale=TRUE, na.action = "na.omit")
      }else{
        pcoa.res=capscale(data.frame(t(matched))~1,distance=dis, na.action = "na.omit")
      }
      # fill background sample positions from PCoA without extra sample
      if(ord.index==0){
        sample.coords[1:lastIndexRef,]=pcoa.res$CA$u[, dimensions]
      }else{
        print(paste("Completed ordination",ord.index))
        # extract coordinates of sample ordinated together with background
        sample.coords[(lastIndexRef+ord.index),]=pcoa.res$CA$u[(lastIndexRef+1), dimensions]
      }
    }
    # check
    #print(sample.coords[1:100,])

  }else{
    # carry out PCoA
    if(dis=="cor"){
      # carry out standard PCA
      # from vegan's doc: using correlation coefficients instead of covariances will give a more balanced ordination (scale=TRUE)
      print("Carrying out PCA...")
      pcoa.res=rda(data.frame(t(abundances)), scale=TRUE, na.action = "na.omit")
    }else{
      pcoa.res=capscale(data.frame(t(abundances))~1,distance=dis, na.action = "na.omit")
    }
    if(!is.null(metadata) && topMetadata>0){

      # carry out envfit
      # Lisa: >> betadisper computes eigenvectors that differ slightly from those in capscale
      # to be consistent, give the vector slot of the output of eigen to envfit
      # from betadisper code:
      # n <- attr(dis, "Size")
      # x <- matrix(0, ncol = n, nrow = n)
      # x[row(x) > col(x)] <- dis^2x <- x + t(x)
      # x <- dblcen(x)
      # e <- eigen(-x/2, symmetric = TRUE) <<
      # for the moment ignored, since betadisper/adonis is only used for assessment of cluster variability & quality and not for plotting
      ef=envfit(pcoa.res,metadata,perm=env.permut, choices=dimensions)
      #print(ef$vectors)
      #print(names(ef$vectors$r))
      if(length(metadata.to.plot)>0){
        for(metadatum.to.plot in metadata.to.plot){
          index.metadatum.to.plot=which(names(ef$vectors$r)==metadatum.to.plot)
          print(paste(sep="","R2 of ",metadatum.to.plot,": ",ef$vectors$r[index.metadatum.to.plot]))
          print(paste(sep="","P-value (not corrected) of ",metadatum.to.plot,": ",ef$vectors$pvals[index.metadatum.to.plot]))
        }
      }
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
  redun.method="PCoA"
  if(dis=="cor"){
    redun.method="PCA"
  }

  if(!addToRefStepwise){
    # proportion of variance explained
    eig.sum=summary(eigenvals(pcoa.res))
    var.explained.1=eig.sum[2,dimensions[[1]]] # second row of eigenvalue summary: proportion of variance explained
    var.explained.2=eig.sum[2,dimensions[[2]]]
    xlab=paste(redun.method,dimensions[1]," [",round(var.explained.1,2),"]",sep="")
    ylab=paste(redun.method,dimensions[2]," [",round(var.explained.2,2),"]",sep="")
  }else{
    xlab=paste(redun.method,dimensions[1],sep="")
    ylab=paste(redun.method,dimensions[2],sep="")
  }

  if(length(colors)>0){
    # use pre-assigned colors
    if(length(colors)!=ncol(abundances)){
      stop("There should be as many colors as samples in the color vector!")
    }
  }else{
    colors=rep("gray",ncol(abundances))
    if(groupAttrib!=""){
      groups=metadata[[groupAttrib]]
      colors=assignColorsToGroups(groups,refName = refName, myColors = groupColors)
    }else if(length(groups)>0){
      colors=assignColorsToGroups(groups, refName = refName, myColors = groupColors)
    }
  }

  # select indices of samples to show (if hiddenSamples is empty, this will be all)
  selected.sample.indices=setdiff(1:ncol(abundances),hiddenSamples)

  # assign cluster memberships as shapes
  if(length(clusters)>0){
    clusnum=length(unique(clusters))
    if(clusnum>length(clus.pch.values)){
      stop(paste("No more than",length(clus.pch.values),"clusters can be visualized."))
    }
    pch.value=c()
    clustersymbol.lookup=list()
    clus.values.counter=1
    for(clus.index in 1:length(clusters)){
      current.clus=clusters[clus.index]
      if(!(current.clus %in% names(clustersymbol.lookup))){
        clustersymbol.lookup[[current.clus]]=clus.pch.values[clus.values.counter]
        clus.values.counter=clus.values.counter+1
      }
    }
    for(clus.index in 1:length(clusters)){
      pch.value=c(pch.value,clustersymbol.lookup[[clusters[clus.index]]])
    }
    pch.value=pch.value[selected.sample.indices]
    #print(pch.value[20:50])
  }

  #print(selected.sample.indices)
  #print(colors)
  # determine range
  if(is.null(xlim)){
    xlim=range(pcoa.res$CA$u[selected.sample.indices,dimensions[1]])
  }
  if(is.null(ylim)){
    ylim=range(pcoa.res$CA$u[selected.sample.indices,dimensions[2]])
  }

  if(!is.null(sizes) || length(sizes)>0){
    display.size.legend=TRUE
    # shift into positive range
    if(min(sizes)<0){
      sizes=sizes-min(sizes)
    }
    # scale between 0.5 and 2.5
    sizes=(sizes-min(sizes))/(0.5*max(sizes)-min(sizes))
    sizes=sizes+min.cex # make sure there is no dot of size zero
    #print(range(sizes))
  }else{
    sizes=rep(par()$cex,length(selected.sample.indices)) # default size
  }

  if(!addToRefStepwise){
    plot(pcoa.res$CA$u[selected.sample.indices,dimensions], cex=sizes[selected.sample.indices], col=colors[selected.sample.indices], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, bg=colors[selected.sample.indices], pch=pch.value, ...)
  }else{
    plot(sample.coords[selected.sample.indices,dimensions], cex=sizes[selected.sample.indices], col=colors[selected.sample.indices], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, bg=colors[selected.sample.indices], pch=pch.value, ...)
  }

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

  groups.factor=factor(groups,levels=unique(groups))
  # add ellipses if requested
  if(drawEllipse==TRUE && !addToRefStepwise){
    if(length(groups)==0){
      warning("Please specify groups in order to draw ellipses!")
    }else{
      # draw transparent ellipses around samples of the same groups with given confidence
      # by default, color order is not correct, because ordiellipse calls factor on the groups, which sorts entries alphabetically
      # for this reason, factor with correct ordering is reassigned to groups
      # color ellipse
      #ordiellipse(pcoa.res,scaling=0,groups=groups.factor,draw="polygon",col=unique(colors),alpha=0.25,conf=ellipseConf,lwd=0.5, kind="sd", border=0)
      # color border rather than ellipse itself
      ordiellipse(pcoa.res,scaling=0,groups=groups.factor,draw="polygon",alpha=0.25,conf=ellipseConf,lwd=0.5, kind="sd", border=unique(colors))
      # adonis fails for step-wise addition of samples to reference
      adonis_results = adonis(data.frame(t(abundances)) ~ groups.factor, permutations = env.permut, method=dis)
      # variance explained through groups
      adonis.r2=adonis_results$aov.tab[1,5]
      adonis.pval=adonis_results$aov.tab[1,6]
      print("Adonis to test for significant difference in group compositions")
      print(paste("Adonis R2: ",round(adonis.r2,4),", p-value: ",round(adonis.pval,4),sep=""))
    }
  }

  # report cluster quality of groups
  if(length(groups)>0){
    if(clusterQualityIndex != "" && clusterQualityIndex != "none"){
      if(clusterQualityIndex=="CH"){
        clusterQualityIndex="Calinski_Harabasz"
      }
      print(paste("Cluster quality index", clusterQualityIndex))
      #print(is.numeric(abundances))
      if(clusterQualityIndex=="silhouette"){
        print(silhouette(abundances,groups=groups,method=dis))
      }else{
        # note that clusterCrit gives an error for large input matrices, but not for small ones; this error is not yet fixed
        print(intCriteria(t(abundances), part=groupsToNumeric(groups), crit=clusterQualityIndex)[[1]])
      }
    }
  }

  if(topTaxa>0){
    # taken from biplot.pcoa in ape
    # columns are taxa
    Y=t(abundances)
    n <- nrow(Y) # n = sample number
    # standardize eigen vectors (subtract mean and divide by standard deviation)
    # eigen vector dimensions: samples x selected dimensions
    ev.stand <- scale(pcoa.res$CA$u[,dimensions])
    # covariance between taxa (as columns) and selected principal components (standardized eigen vectors)
    S <- cov(Y, ev.stand)
    # scale S by the eigen values
    U <- S %*% diag((pcoa.res$CA$eig[dimensions]/(n-1))^(-0.5))
    colnames(U) <- colnames(pcoa.res$CA$u[,dimensions])
    # U dimensions: taxa x selected dimensions

    # select top covarying taxa in U
    # arrrow length codes strength of covariance with eigen vectors
    norms=apply(U,1,myNorm)
    sorted=sort(norms,index.return=TRUE,decreasing=TRUE)
    sorted.top.indices=sorted$ix[1:topTaxa]
    pvalues=c()
    U.sub=U[sorted.top.indices,]
    Y.sub=Y[,sorted.top.indices]
    was.permuted=FALSE
    # if requested, carry out permutations
    if(!is.null(permut) && !is.na(permut) && permut>0){
      was.permuted=TRUE
      # taxa vs permutations
      permuted.norms=matrix(NA,nrow=ncol(Y.sub),ncol=permut)
      # carry out permutation test
      for(iteration in 1:permut){
        # shuffle abundances separately per column (in Y, taxa are columns)
        Y.rand=apply(Y.sub,2,base::sample)
        # only look at top covarying taxa
        S.rand=cov(Y.rand,ev.stand)
        # scale the shuffled S by eigen values
        U.rand <- S.rand %*% diag((pcoa.res$CA$eig[dimensions]/(n-1))^(-0.5))
        # taxon arrows are defined by U, which has as many rows as taxa and as many columns as selected eigen vectors
        # compute arrow length as vector norm for each row (each row represents one arrow)
        permuted.norms[,iteration]=apply(U.rand,1,myNorm)
      }
      # compute signficance of vector norms
      for(taxon.index in 1:ncol(Y.sub)){
        obs.norm=norms[taxon.index]
        # check how many permuted norms are larger than observed norm
        r=length(which(permuted.norms[taxon.index,]>obs.norm))
        # compute parameter-free p-value
        pvalues=c(pvalues,(r+1)/(permut+1))
      }
      # adjust p-values for multiple testing and discard corrected p-values below selected significance level
      pvalues=p.adjust(pvalues,method=pAdjMethod)
      sig.pvalue.indices=which(pvalues<qvalThreshold)
      # only keep significant top covarying taxa
      U.selected=U.sub[sig.pvalue.indices,]
    } # end permutation test
    else{
      U.selected=U.sub
      sig.pvalue.indices=c()
    }
    if(length(sig.pvalue.indices)>0 || !was.permuted){
      if(was.permuted){
        print(paste("Among the top ",topTaxa," covarying taxa, ",length(sig.pvalue.indices)," are significant.",sep=""))
        for(sig.taxon.name in rownames(U.selected)){
          print(sig.taxon.name)
        }
      }
      if(length(hiddenTaxa)>0){
        visibleTaxa.indices=c()
        # loop selected taxa
        for(U.taxon.index in 1:nrow(U.selected)){
          if(!(rownames(U.selected)[U.taxon.index] %in% hiddenTaxa)){
            visibleTaxa.indices=c(visibleTaxa.indices,U.taxon.index)
          } # check whether selected taxon should be hidden
        }
        U.selected = U.selected[visibleTaxa.indices,]
      }
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
    }else{
      if(was.permuted){
        warning("None of the top-covarying taxa is significant.")
      }
    }
  }

  if(!is.null(metadata) && topMetadata>0){
    if(length(indices.vectors)>0 || length(metadata.to.plot)>0){
      if(length(metadata.to.plot) > 0){
        #print(names(ef$vectors$r[pvals.vectors.sorted$ix]))
        for(metadatum.to.plot in metadata.to.plot){
          index.metadatum.to.plot=which(names(ef$vectors$r[pvals.vectors.sorted$ix])==metadatum.to.plot)
          indices.vectors=c(indices.vectors,index.metadatum.to.plot)
        }
      }
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
      # significance of centroid separation: now assessed with adonis
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
    # cex=0.9
    legend("topright",legend=unique(groups[selected.sample.indices]),cex=0.6, pch = rep("*",length(unique(groups[selected.sample.indices]))), col = unique(colors[selected.sample.indices]), bg = "white", text.col="black")
  }
  if(length(clusters)>0){
    legend("topleft",legend=unique(clusters[selected.sample.indices]),cex=0.9, pch = unique(pch.value), col = "black", bg = "white", text.col="black")
  }
  # 1 = default size
  if(display.size.legend==TRUE){
    min.legend="min"
    max.legend="max"
    if(size.legend!=""){
      min.legend=paste(min.legend," ",size.legend,sep="")
      max.legend=paste(max.legend," ",size.legend,sep="")
    }
    legend("bottomleft",legend=c(min.legend,max.legend),pt.cex=c(min.cex,max.cex), pch = c(1,1), col = "black", bg = "white", text.col="black")
  }

  # A non-significant result in betadisper is not necessarily related to a significant/non-significant result in adonis.
  # Betadisper tests homogeneity of dispersion among groups, which is a condition (assumption) for adonis.
  # Betadisper can be done to see if one group has more compositional variance than another.
  # Adonis tests whether composition among groups is similar or not.
  if(groupDispersion){
    beta.out=betadisper(vegdist(t(abundances),method=dis), groups.factor, type='centroid')
    print("Tukey's HSD to test for significant differences in group variance (betadisper)")
    t.hsd=TukeyHSD(beta.out)
    print(t.hsd$group)
    # fourth column: p-values
    if(length(which(t.hsd$group[,4]<0.05))){
      print("Variance differs signficantly for at least 2 groups, so adonis results may be biased!")
    }
    boxplot(beta.out,xlab="Group")
  }

}

#' @title Assign colors to groups
#'
#' @description Given a group membership vector, each group receives
#' a unique color such that the resulting color vector is as long as
#' the group membership vector.
#'
#' @param groups a vector of group memberships
#' @param refName the name of the reference group; receives gray as a color
#' @param myColors a map of predefined colors for groups
#' @param returnMap whether to return the color map together with the colors
#' @return a vector of colors; if returnMap is true, a vector of colors and the color map
#' @export
#'
assignColorsToGroups<-function(groups, refName="ref", myColors = NULL, returnMap=FALSE){
  refContained=FALSE
  if(refName %in% groups){
    refContained=TRUE
  }
  groupNum=length(unique(groups))
  if(refContained){
    groupNum=groupNum-1 # reference counts extra
  }
  col.vec = seq(0,1,1/groupNum)
  hues = hsv(col.vec)
  #print(hues)
  hueCounter=1
  colors=c()
  # fill the color map
  if(is.null(myColors)){
    myColors=list()
    for(group.index in 1:length(groups)){
      group=as.character(groups[group.index])
      if(!(group %in% names(myColors))){
        if(group==refName){
          myColors[[group]]="gray"
        }else{
          myColors[[group]]=hues[hueCounter]
          hueCounter=hueCounter+1
        }
      }
    } # loop samples
  } # color map already filled

  #print(myColors)
  # assign colors from the color map
  for(group.index in 1:length(groups)){
    #print(groups[group.index])
    #print(myColors[[as.character(groups[group.index])]])
    colors=c(colors,myColors[[as.character(groups[group.index])]])
  }

  if(returnMap){
    res=list(colors,myColors)
    names(res)=c("colors","colorMap")
    return(res)
  }
  return(colors)
}

# convert a group vector with strings in a numeric vector of cluster membership integers
groupsToNumeric<-function(groups=c()){
  clus.mem=groups
  groups.unique=unique(groups)
  for(i in 1:length(groups.unique)){
    current.group=groups.unique[i]
    group.indices=which(groups==current.group)
    clus.mem[group.indices]=i
  }
  return(as.integer(as.numeric(clus.mem)))
}


# compute the norm of a vector
myNorm<-function(x){
  if(!is.numeric(x)){
    stop("x should be a numeric vector.")
  }
  return(sqrt(sum(x^2)))
}

# Rarefaction combined with sample filtering
#
# Rarefy a matrix to the given minimum count number column-wise
# using vegan's rrarefy function. If columns have less than the minimum count number,
# they are discarded. Rows that have a sum of zero after rarefaction are also discarded.
# x a matrix
# min minimum count to which x is to be rarefied (if equal to zero, the minimum column sum is taken as min)
# a list with the rarefied matrix (rar) and the indices of the columns that were kept (colindices)
rarefyFilter<-function(x,min = 0){
  keep=c()
  if(min < 0){
    stop("Min should be either 0 or positive.")
  }
  if(min == 0){
    min=min(colsums=apply(x,2,sum))
    print(paste("Rarefy to minimum count",min))
    keep=c(1:ncol(x))
  }else{
    colsums=apply(x,2,sum)
    # there are columns below the minimum
    if(min(colsums) < min){
      # loop column sums
      for(j in 1:ncol(x)){
        if(colsums[j] >= min){
          keep=c(keep,j)
        }
      }
      print(paste("Number of columns",ncol(x)))
      print(paste("Keeping ",length(keep)," columns with column sums equal or above",min))
      x=x[,keep]
    }
  }
  rar=t(vegan::rrarefy(t(x),min))
  zero.indices=which(rowSums(rar)==0)
  # discard taxa with zero sums
  if(length(zero.indices)>0){
    keep.nonzero=setdiff(1:nrow(rar),zero.indices)
    rar=rar[keep.nonzero,]
  }
  res=list(rar,keep)
  names(res)=c("rar","colindices")
  return(res)
}

