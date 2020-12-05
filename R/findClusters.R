#' @title Cluster taxon- or sample-wise
#'
#' @description Cluster sequencing data taxon- or sample-wise.
#' By default, data are clustered sample-wise. The default method is
#' Dirichlet-Multinomial mixtures (DMM) using the DirichletMultinomial package.
#' Note that DMM expects counts and fails if there are taxa that are absent across all samples.
#' Consequently, abundances are scaled to counts and taxa absent across all samples are removed.
#' Also note that counts should not be too large for DMM (not above 10000). If this is the case,
#' they will be scaled (divided by a constant and rounded). In addition, for DMM, the total number
#' of counts should be the same across samples.
#' Partitioning around medoids (PAM) does not expect counts, can deal with absent taxa and large counts.
#' The cluster number can be omitted, in which case the best cluster number is determined in the range
#' provided using either silhouette or a cluster quality index implemented in clusterCrit's function intCriteria.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param method clustering method, supported are dmm and pam
#' @param k cluster number, can be set to NA if cluster number is to be determined with a quality index
#' @param minK the minimum cluster number to test
#' @param maxK the maximum cluster number to test
#' @param qualityIndex by default silhouette computed with silhouette function on bray curtis, can also be CH, Dunn, and other quality indices supported by clusterCrit's intCriteria function or none
#' @return a cluster membership vector
#' @examples
#' data("ibd_taxa")
#' # unequal sample sums in the example data set are due to the removal of low-prevalence taxa
#' clusters=findClusters(ibd_taxa, k=NA, method="pam", qualityIndex="CH")
#' table(clusters)
#' @export
#'

findClusters<-function(abundances, method="dmm", k=3, minK=2, maxK=10, qualityIndex="silhouette"){
  # DMM needs not too large counts and doesn't work when taxa are zero everywhere.
  if(method=="dmm"){
    taxonSums=rowSums(abundances)
    absent.taxa=which(taxonSums==0)
    colSums=colSums(abundances)
    if(length(unique(colSums))>1){
      warning("The total counts across samples are not the same. DMM may be biased.")
    }
    if(length(absent.taxa)>0){
      warning("DMM does not accept taxa with no occurrence across samples. These have been removed.")
      taxa.to.keep=setdiff(1:nrow(abundances),absent.taxa)
      abundances=abundances[taxa.to.keep,]
    }
    if(!is.Count.Matrix(abundances)){
      warning("DMM only accepts counts, thus abundances have been scaled to counts.")
      vec=as.vector(abundances)
      scaling.factor=1/min(vec[vec!=0])
      if(scaling.factor>10000){
        scaling.factor=10000
      }
      print(paste("Abundances have been scaled with factor:",scaling.factor))
      abundances=round(abundances*scaling.factor)
    }else{
      if(max(abundances)>10000){
        scaling.factor=max(abundances)/10000
        warning(paste("Counts are too large for DMM and are scaled by ",scaling.factor,", followed by rounding.",sep=""))
        abundances=round(abundances/scaling.factor)
      }
    }
  }

  if(is.na(k) && (qualityIndex=="none" || qualityIndex=="")){
    stop("When k is not given, a cluster quality index has to be provided.")
  }

  if(qualityIndex=="CH"){
    qualityIndex="Calinski_Harabasz"
  }

  # k provided
  if(!is.na(k) && !is.null(k)){
    groups=findClustersGivenK(abundances, method=method,k=k)
    if(qualityIndex!="none"){
      clusterQuality=NA
      if(qualityIndex=="silhouette"){
        clusterQuality=silhouette(abundances,groups=groups,method="bray")
      }else{
        # note that clusterCrit gives an error for large input matrices, but not for small ones; this error is not yet fixed
        clusterQuality=intCriteria(t(abundances),groups,c(qualityIndex))
      }
      print(paste("Cluster quality according to",qualityIndex,":",clusterQuality))
    }
    return(groups)
  }else{
    if(qualityIndex=="none"){
      stop("Please provide a quality index.")
    }
    # store group memberships
    kToTest=c(minK : maxK)
    groupMatrix=matrix(nrow=length(kToTest),ncol=ncol(abundances))
    qualityValues=c()
    rowCounter=1
    for(k in kToTest){
      groupMatrix[rowCounter,]=findClustersGivenK(abundances, method=method,k=k)
      qualityValue=NA
      if(qualityIndex=="silhouette"){
        qualityValue=silhouette(abundances,groups=groupMatrix[rowCounter,],method="bray")
      }else{
        qualityValue=intCriteria(t(abundances),groupMatrix[rowCounter,],c(qualityIndex))
      }
      print(paste("k equal",k,"has quality value",qualityValue))
      qualityValues=c(qualityValues,qualityValue)
      rowCounter=rowCounter+1
    }
    #print(qualityValues)
    unlisted=unlist(qualityValues)
    na.indices=which(is.na(unlisted))
    # deal with missing values
    if(length(na.indices)>0){
      non.na=setdiff(c(1:length(unlisted)),na.indices)
      unlisted=unlisted[non.na]
      kToTest=kToTest[non.na]
    }
    index.best=bestCriterion(unlisted, crit=qualityIndex) # also works for silhouette since the name is the same
    #print(index.best)
    print(paste("The optimal cluster number according to",qualityIndex,"is",kToTest[index.best],"with value",qualityValues[index.best]))
    return(as.vector(groupMatrix[index.best,]))
  }
}

# find clusters given a particular value for k
findClustersGivenK<-function(abundances, method="dmm", k=5){
  groups=c()
  if(method=="dmm"){
    dmn.out=dmn(t(abundances),k=k, verbose=TRUE)
    # helps cluster number selection
    # BIC/AIC: the smaller, the better
    print("DMM goodness of fit:")
    print(goodnessOfFit(dmn.out))
    clus.assignment=mixture(dmn.out)
    # assign DMM type
    for(row.index in 1:nrow(clus.assignment)){
      clus.contribs=clus.assignment[row.index,]
      dmm.type=which(clus.contribs==max(clus.contribs))
      groups=c(groups,dmm.type)
    }
  }else if(method=="pam"){
    pam.out=pam(t(abundances),k=k)
    groups=pam.out$clustering
  }
  return(groups)
}


#' @title Compute transition probabilities between clusters
#'
#' @description The cluster vector assigns to each
#' sample its cluster membership. If
#' samples are in temporal order,
#' transition probabilities between
#' clusters can be computed. A group membership
#' vector can be optionally provided to prevent
#' computation of transitions between different
#' groups. This is useful when time series were
#' collected for several related experiments
#' (e.g. study participants). The order of samples
#' in the cluster and group membership vectors is
#' supposed to be the same.
#' If logodds is true, the log ratio between the observed and
#' expected transition probability is computed, where the latter
#' is the product of the prior probabilities for the two clusters. A log
#' ratio above 1 or below -1 means that observed transition probabilities are
#' greater or smaller than expected at random.
#' If a binary metadata item is provided, the function
#' counts how often transitions occur within a sliding window
#' with and without a metadata change. In this case, two matrices are
#' returned: one with transition frequencies in the presence of
#' a metadata change and one with transition frequencies in the
#' absence of that change. In this case, intra-cluster
#' transitions are not counted.
#' Note that within the window, temporal order between
#' metadata and cluster change is not enforced (so the
#' cluster may change before the metadata changes).
#'
#' @param clus.vec the vector with cluster memberships
#' @param groups an optional vector with group memberships
#' @param metadata.vec an optional binary metadata item
#' @param windowSize in case metadata are provided: the size of the sliding window
#' @param freq if true, transition frequencies instead of probabilities are returned (for metadata.vec, always true)
#' @param logodds if true, compute logarithm of the ratio of the observed and the expected transition probability
#' @return a matrix or, if metadata.vec provided, a list with 2 matrices, the first without metadata change and the second with metadata change
#' @examples
#' # generate random cluster memberships
#' clus.vec=round(runif(100,min=1,max=5))
#' trans.mat=transitionProbabs(clus.vec)
#' # display transition matrix as a network with igraph
#' gr=graph_from_adjacency_matrix(trans.mat,mode="directed",weighted = TRUE)
#' # plot graph with transition probabilities as edge thickness
#' edge.weights=as.vector(t(trans.mat))
#' edge.weights=edge.weights[edge.weights>0]
#' plot(gr,edge.width=edge.weights*10) # 10 is a factor to scale edge weights, may have to be adjusted
#' @export
transitionProbabs<-function(clus.vec=c(), groups=c(), metadata.vec=c(), windowSize=4, freq=FALSE, logodds=FALSE){
  na.clus=which(is.na(clus.vec))
  if(length(na.clus)>0){
    stop("Missing values in the cluster assignment are not allowed.")
  }
  if(logodds && length(metadata.vec)>0){
    warning("Log odds ratio is not calculated when metadata are provided.")
    logodds=FALSE
  }
  unique.groups=unique(clus.vec)
  clus.transitions=matrix(0,nrow=length(unique.groups),ncol=length(unique.groups))
  clus.triggered.transitions=matrix(0,nrow=length(unique.groups),ncol=length(unique.groups))
  rownames(clus.transitions)=unique.groups
  colnames(clus.transitions)=rownames(clus.transitions)
  rownames(clus.triggered.transitions)=unique.groups
  colnames(clus.triggered.transitions)=rownames(clus.transitions)

  # metadata provided
  if(length(metadata.vec)>0){
    for(i in 1:(length(clus.vec)-windowSize)){
      forbidden=FALSE
      window.clusvec=clus.vec[i:(i+windowSize)]
      window.metadata.vec=metadata.vec[i:(i+windowSize)]
      # omit windows spanning two groups
      if(length(groups)>0){
        window.groups=groups[i:(i+windowSize)]
        if(length(unique(window.groups))>1){
          forbidden=TRUE
        }
      }
      if(!forbidden){
        clus.values=unique(window.clusvec)
        metadata.values=unique(window.metadata.vec)
        metadata.values=setdiff(metadata.values,NA) # remove NA
        clus.change=FALSE
        metadata.change=FALSE
        # cluster transition
        if(length(clus.values)>1){
          clus.change=TRUE
        }
        # metadata change
        if(length(metadata.values)>1){
          metadata.change=TRUE
        }
        if(clus.change){
          clus.index1=which(rownames(clus.transitions)==clus.values[1])
          clus.index2=which(rownames(clus.transitions)==clus.values[2])
          if(metadata.change){
            clus.triggered.transitions[clus.index1,clus.index2]=clus.triggered.transitions[clus.index1,clus.index2]+1
          }else{
            clus.transitions[clus.index1,clus.index2]=clus.transitions[clus.index1,clus.index2]+1
          }
        } # cluster membership change
      } # not forbidden
    } # loop over windows
  }else{
    prev.clus=clus.vec[1]
    prev.group=NA
    if(length(groups)>0){
      prev.group=groups[1]
    }
    transitions=0
    forbidden=FALSE
    for(i in 2:length(clus.vec)){
      forbidden=FALSE
      clus=clus.vec[i]
      if(length(groups)>0){
        group=groups[i]
        # only count transitions within 1 group
        if(group!=prev.group){
          forbidden=TRUE
        }
        prev.group=group
      }
      if(!forbidden){
        index.prevclus=which(rownames(clus.transitions)==prev.clus)
        index.clus=which(rownames(clus.transitions)==clus)
        transitions=transitions+1
        clus.transitions[index.prevclus,index.clus]=clus.transitions[index.prevclus,index.clus]+1
      }
      prev.clus=clus
    } # end loop clus.vec

    # convert into probabilities when freq is false or logodds is true
    if(!freq || logodds){
      clus.transitions=clus.transitions/transitions
    }
    # compute log-ratio of observed probability and prior probability
    # prior probability is the number of times a cluster is observed divided by sample number
    if(logodds){
      prior.probabs=table(clus.vec)/length(clus.vec)
      print(paste("Prior probabilities"))
      print(prior.probabs)
      prior.probabs.ordered=prior.probabs
      # give prior probabs the same order as in transition matrix
      group.counter=1
      for(group in unique.groups){
        index.prior.probabs=which(rownames(prior.probabs)==group)
        prior.probabs.ordered[group.counter]=prior.probabs[index.prior.probabs]
        group.counter=1+group.counter
      }
      # compute expected transition probabilities
      expected.clus.transitions = matrix(0,nrow=length(unique.groups),ncol=length(unique.groups))
      for(i in 1:length(unique.groups)){
        for(j in 1:length(unique.groups)){
          expected.clus.transitions[i,j]=prior.probabs.ordered[i]*prior.probabs.ordered[j]
        }
      }
      # compute logodds matrix
      clus.transitions=log(clus.transitions/expected.clus.transitions)
      clus.transitions[is.infinite(clus.transitions)]=NA
    }
  }
  if(length(metadata.vec)>0){
    res=list(clus.transitions,clus.triggered.transitions)
    names(res)=c("transit","trigger.transit")
    return(res)
  }
  return(clus.transitions)
}

# check whether a matrix has only integers.
# missing values are not treated.
is.Count.Matrix<-function(abundances){
  res=sapply(abundances,is.Count)
  if(length(which(res==FALSE))>0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

# check whether a number is an integer
is.Count<-function(x){
  # integer modulo 1 should always be 0
  if(x%%1==0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
