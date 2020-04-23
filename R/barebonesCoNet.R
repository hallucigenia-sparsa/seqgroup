#' @title Barebones R implementation of CoNet
#'
#' @description Build a network using Pearson, Spearman, Kullback-Leibler and/or Bray-Curtis.
#' The original CoNet implementation with extended functionality is available at: \href{http://psbweb05.psb.ugent.be/conet/}{http://systemsbiology.vub.ac.be/conet}.
#' For export of data to the original CoNet, see \code{\link{exportToCoNet}}, for a generic network building function wrapping
#' barebonesCoNet and other network inference methods, see \code{\link{buildNetwork}}.
#'
#' @details If renorm and permutandboot are both set to TRUE, p-value computation is equal to the ReBoot procedure implemented
#' in CoNet. If more than one method is selected and p-value computation is enabled, p-values are merged with Fisher's method, multiple
#' testing correction (if enabled) is applied on the merged p-value and the merged p-value is reported. If p-value computation is not enabled,
#' the method number is reported as association strength. Note that for a single dissimilarity method, weights are scaled to the absolute distance
#' from the mean value (for bounded Bray Curtis, this is 0.5 and for KLD, this is the mean of the observed scores), such that a larger edge weight
#' means a stronger association.
#' Edge signs (co-presence/mutual exclusion) are assigned using thresholds (T.up/T.down directly or indirectly via top/bottom initial edge number).
#' Co-presence (high correlation/low dissimilarity) is encoded in green, mutual exclusion (low correlation/high dissimilarity) in red and sign conflicts
#' (lack of agreement between methods) in gray.
#' When metadata are provided, a bipartite network is computed. To circumvent bipartite network computation, metadata can be appended to abundances to form a single
#' input matrix, but in this case, preprocessing on abundance data needs to be carried out before.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param metadata an optional data frame with metadata items as columns, where samples are in the same order as in abundances and all items are numeric; a bipartite network will be computed
#' @param methods network construction methods, values can be combinations of: "pearson", "spearman", "kld" or "bray"
#' @param T.up upper threshold for scores (when more than one network construction method is provided, init.edge.num is given and/or p-values are computed, T.up is ignored)
#' @param T.down lower threshold for scores (when more than one network construction method is provided, init.edge.num is given and/or p-values are computed, T.down is ignored)
#' @param method.num.T threshold on method number (only used when more than one method is provided)
#' @param pval.T threshold on p-value (only used when permut, permutandboot or pval.cor is true); if several methods are provided, only applied after merge
#' @param init.edge.num the number of top and bottom initial edges (init.edge.num overrides T.up/T.down, set to NA to respect T.up/T.down for a single method)
#' @param min.occ only keep rows with at least the given number of non-zero values (carried out before network construction)
#' @param keep.filtered sum all filtered rows and add the sum vector as additional row
#' @param norm normalize matrix (carrried out after filtering)
#' @param stand.rows standardize rows by dividing each entry by its corresponding row sum, applied after normalization
#' @param pval.cor compute p-values of correlations with cor.test (only valid for correlations; takes precedence over permut and permutandboot with or without renorm)
#' @param permut compute p-values on edges with a permutation test
#' @param renorm use renormalization when computing permutation distribution (only applied to correlations; cannot be combined with metadata)
#' @param permutandboot compute p-values from both permutation (with or without renorm) and bootstrap distribution
#' @param iters number of iterations for the permutation and bootstrap distributions
#' @param bh multiple-test-correct using Benjamini-Hochberg; if several methods are provided bh is applied to merged p-value
#' @param pseudocount count added to zeros prior to taking logarithm (for KLD, p-value merge and significance)
#' @param plot plot score or, if permut, permutandboot or pval.cor is true, p-value distribution, in both cases after thresholding
#' @param verbose print the number of positive and negative edges and, if permut, permutandboot or pval.cor is true, details of p-value computation
#' @return igraph object with absolute association strengths, number of supporting methods or, if permut, permutandboot or pval.cor is true, significances (-1*log10(pval)) as edge weights
#' @examples
#' data("ibd_taxa")
#' data("ibd_lineages")
#' ibd_genera=aggregateTaxa(ibd_taxa,lineages = ibd_lineages,taxon.level = "genus")
#' min.occ=nrow(ibd_genera)/3
#' # p-values for the 50 strongest positive and 50 strongest negative Spearman correlations
#' cn=barebonesCoNet(ibd_genera,methods="spearman",init.edge.num=50,min.occ=min.occ,permutandboot=TRUE)
#' plot(cn)
#' # combine Bray Curtis and Spearman and threshold on method number
#' plot(barebonesCoNet(ibd_genera,methods=c("spearman","bray"),init.edge.num=50,min.occ = min.occ))
#' @export
barebonesCoNet<-function(abundances, metadata=NULL, methods=c("spearman","kld"), T.up=NA, T.down=NA, method.num.T=2, pval.T=0.05, init.edge.num=max(2,round(sqrt(nrow(abundances)))), min.occ=0, keep.filtered=TRUE, norm=FALSE, stand.rows=FALSE, pval.cor=FALSE, permut=FALSE, renorm=FALSE, permutandboot=FALSE, iters=100, bh=TRUE, pseudocount=0.00000000001, plot=FALSE, verbose=FALSE){

  correlations=c("pearson","spearman")
  bh.selected=bh
  pval.T.selected=pval.T
  N=nrow(abundances)
  forbidden.combis=NULL
  sign.matrix.list=list()
  resultList=list()
  res.graph=NULL
  total.edge.num=(N*(N-1))/2
  #print(total.edge.num)
  #print(init.edge.num)
  default.init.edge.num=max(2,round(sqrt(N)))
  default.pseudocount=0.00000000001 # used to convert p-values to significances to avoid loosing zero p-values (cannot be modified by the user)
  taxon.names=c()
  metadata.names=c()

  methods=unique(methods)

  print(paste("Network construction with method(s):",paste0(methods,collapse=", ")))

  ### Checks
  if(length(methods)>1){
    print("Upper and lower threshold are ignored when more than one method is selected.")
    T.up=NA
    T.down=NA
    bh=FALSE # multiple-testing correction is applied after the merge
    pval.T=1 # p-values are filtered after merge
    if(is.na(init.edge.num) || init.edge.num<1){
      warning(paste("When more than 1 method is selected, init.edge.num has to be provided! init.edge.num is set now to 1/3 of the total edge number, namely",default.init.edge.num))
      init.edge.num=default.init.edge.num
    }
  }

  if(!is.null(metadata) && renorm==TRUE){
    stop("Renormalisation with metadata is not supported.")
  }

  if((permut==TRUE || permutandboot==TRUE) && iters < 1){
    stop("iters should be at least 1.")
  }

  if(!is.na(T.up) || !is.na(T.down)){
    init.edge.num=NA
  }

  if(!is.na(init.edge.num) && total.edge.num < init.edge.num){
    stop(paste("The total edge number of ",total.edge.num," is smaller than the initial edge number. Please reduce the initial edge number."))
  }
  if(length(methods)==1){
    method=methods[1]
    if(pval.cor == TRUE & (method == "kld" || method == "bray")){
      stop("P-value computation with cor.test is only possible for correlations!")
    }
  }

  if(pseudocount<=0){
    stop("Pseudocount has to be set to a small positive number.")
    #pseudocount=min(matrix[matrix>0])/100
    #print(paste("Pseudocount:",pseudocount))
  }

  ### Preprocessing
  abundances = filterTaxonMatrix(abundances,keepSum = keep.filtered, minocc=min.occ)

  # normalize matrix
  if(norm == TRUE){
    abundances=normalize(abundances)
  }
  # normalize row-wise
  if(stand.rows == TRUE){
    abundances = t(normalize(t(abundances)))
  }

  if(!is.null(metadata)){
    metadata.names=colnames(metadata) # metadata items are columns
    taxon.names=rownames(abundances)
    #print(dim(metadata))
    abundances=rbind(abundances,t(metadata)) # append metadata
  }

  # update row number
  N=nrow(abundances)

  ### Network construction
  # loop selected methods
  for(method in methods){
    sign.matrix=matrix(NA,nrow=N,ncol=N)
    renorm.selected=renorm
    print(paste("Processing method",method))
    # initial edge number provided
    if(!is.na(init.edge.num) && init.edge.num>0){
      #  get score matrix
      scores=getScores(abundances,method=method, pseudocount=pseudocount)
      #print(length(as.vector(scores[lower.tri(scores)])))
      # sort score vector representing lower triangle of score matrix in ascending order
      # note that lower.tri takes values column-wise, not row-wise
      scorevec=sort(as.vector(scores[lower.tri(scores)]),decreasing = FALSE)
      #print(scorevec[1:10])
      # determine lower threshold
      T.down=scorevec[init.edge.num]
      # determine upper threshold
      T.up=scorevec[(length(scorevec)-init.edge.num)]
      print(paste("Lower threshold for initial edge number (",init.edge.num,"): ",T.down,sep=""))
      print(paste("Upper threshold for initial edge number (",init.edge.num,"): ",T.up,sep=""))
      # compute sign matrix
      sign.matrix.list[[method]]=computeSignMatrix(scores,method)
      # set scores to avoid to NA
      forbidden.combis=scores
      forbidden.combis[forbidden.combis>T.up]=Inf
      forbidden.combis[forbidden.combis<T.down]=Inf
      forbidden.combis[!is.infinite(forbidden.combis)]=NA
      #forbidden.indices=which(forbidden.combis>T.down && forbidden.combis<T.up,arr.ind = TRUE)
      scores[is.na(forbidden.combis)]=NA
      forbidden.combis=scores
      if(!is.null(metadata)){
        forbidden.combis=forbiddenCombisBipartite(separator.index = (length(taxon.names)+1), forbidden.combis = forbidden.combis)
      }
      #print(scores)
      # thresholds are not applied when initial edge number is indicated
      T.down=NA
      T.up=NA
      if(renorm == TRUE & (method == "kld" || method == "bray")){
        renorm.selected=FALSE # kld and bray are compositionality-robust
      }
    } # end initial edge number given
    else{
      if(!is.null(metadata)){
        forbidden.combis=getScores(abundances,method=method, pseudocount=pseudocount)
        # compute sign matrix
        sign.matrix.list[[method]]=computeSignMatrix(forbidden.combis,method)
        forbidden.combis=forbiddenCombisBipartite(separator.index = (length(taxon.names)+1), forbidden.combis = forbidden.combis)
      }
    }
    res=computeAssociations(abundances,method=method, forbidden.combis = forbidden.combis, pval.T = pval.T, bh=bh, T.down=T.down, T.up=T.up, pval.cor = pval.cor, renorm=renorm.selected, permut=permut, permutandboot = permutandboot, verbose=verbose, plot=plot, iters=iters, pseudocount=pseudocount)
    resultList[[method]]=res
    # no initial edge number specified: collect sign matrix
    if((is.na(init.edge.num) || init.edge.num==0) && is.null(metadata)){
      sign.matrix.list[[method]]=res$signs
    }
  }

  #print(length(sign.matrix.list))

  if(!is.null(metadata)){
    print(paste("Associations computed for",length(taxon.names),"taxa and",length(metadata.names),"metadata."))
  }else{
    print(paste("Associations computed for",N,"taxa."))
  }

  score.matrix=matrix(0,nrow=N,ncol=N) # stores weights and signs
  pvalue.matrix=matrix(NA,nrow=N,ncol=N) # stores p-values
  sign.matrix=matrix(NA,nrow=N,ncol=N) # stores signs

  # single method
  if(length(methods)==1){
    res=resultList[[methods[1]]]
    sign.matrix=sign.matrix.list[[methods[1]]]
    # make sure nodes have names
    colnames(res$pvalues)=rownames(abundances)
    colnames(res$scores)=rownames(abundances)
    if(permut==TRUE || permutandboot==TRUE || pval.cor==TRUE){
      # for a single method, multiple-testing correction was carried out previously
      # avoid taking the logarithm of zero p-value, else we loose them
      res$pvalues[res$pvalues==0]=default.pseudocount
      pvalue.matrix=res$pvalues
      #print(pvalue.matrix[!is.na(pvalue.matrix)])
      # convert to significances
      sig.matrix=-1*log10(pvalue.matrix)
      # set missing values to absent edges
      sig.matrix[is.na(sig.matrix)]=0
      sig.matrix[is.infinite(sig.matrix)]=0
      diag(sig.matrix)=0 # no self-loops
      #print(sig.matrix[sig.matrix>0])
      # only the lower triangle will be used (this is important, because BH correction is only applied to the lower triangle of the p-value matrix)
      res.graph=graph_from_adjacency_matrix(sig.matrix,mode="lower",weighted=TRUE)
    }else{
      # for a single method, score filtering was carried out previously
      score.matrix=res$scores
      #print(dim(score.matrix))
      # set missing values to absent edges
      score.matrix[is.na(score.matrix)]=0
      #print(score.matrix)
      diag(score.matrix)=0 # no self-loops
      # scale edge weights
      if(method %in% correlations){
        score.matrix=abs(score.matrix)
      }else if(method=="bray"){
        score.matrix=abs(score.matrix-0.5)
      }else if(method=="kld"){
        mean.score=mean(range(score.matrix))
        score.matrix=abs(score.matrix-mean.score)
      }
      res.graph=graph_from_adjacency_matrix(score.matrix,mode="undirected",weighted=TRUE)
    }
  # multiple methods - merge
  }else{
    print("Merging methods...")
    isFirst=TRUE
    #current.sign.matrix=NULL
     for(method in methods){
       res=resultList[[method]]
       # make sure nodes have names
       colnames(res$pvalues)=rownames(abundances)
       colnames(res$scores)=rownames(abundances)
       # avoid taking the logarithm of zero p-value, else we loose them
       res$pvalues[res$pvalues==0]=default.pseudocount
       logpvalue.matrix=log(res$pvalues)
       # after log, set missing values to 0, so they do not set the entire calculation to NA (method is not counted)
       logpvalue.matrix[is.na(logpvalue.matrix)]=0
       # multiply logarithm of individual p-values
       if(permut==TRUE || permutandboot==TRUE || pval.cor==TRUE){
          if(isFirst){
           # print(res$pvalues)
            pvalue.matrix=logpvalue.matrix
          }else{
            pvalue.matrix=pvalue.matrix+logpvalue.matrix
          }
         isFirst=FALSE
       }
       if(isFirst){
         sign.matrix=sign.matrix.list[[method]]
       }else{
         copy.sign.matrix=sign.matrix
         current.sign.matrix=sign.matrix.list[[method]]
         # transfer signs from current matrix
         sign.matrix[current.sign.matrix>0]=1
         sign.matrix[current.sign.matrix<0]=-1
         # check for sign conflicts:  +/+=OK, -/-=OK, -/+=conflict
         # default multiplicator: element-wise matrix multiplication on sign matrix before update
         current.sign.matrix=copy.sign.matrix*current.sign.matrix
         # conflicts: set to NA (will receive a gray color)
         sign.matrix[current.sign.matrix<0]=NA
       }
       # edge-specific method number is needed also for p-value merge (df)
       # add up method number for edges
       res$scores[!is.na(res$scores)]=1
       res$scores[is.na(res$scores)]=0
       #print(length(which(res$scores>0)))
       score.matrix=score.matrix+res$scores
     } # end methods loop

    #print(pvalue.matrix)

    #  filtering on method numbers and graph object creation
    if(!permut && !pval.cor && !permutandboot){
      score.matrix[score.matrix<method.num.T]=0
      diag(score.matrix)=0 # no self-loops
      res.graph=graph_from_adjacency_matrix(score.matrix,mode="undirected",weighted=TRUE)
    }else{
      # correction & filtering on p-values and graph object creation
      #  carry out Fisher's method
      pvalue.matrix=-2*pvalue.matrix
      pvalue.matrix[pvalue.matrix==0]=NA # 0 means missing edge
      # df is the score matrix with the method number per edge
      # higher df will lower the p-value
      pvalue.matrix=pchisq(pvalue.matrix,df=2*score.matrix) # pvalue from chi2 distibution
      #print("Merged p-value matrix")
      #print(pvalue.matrix)
      if(bh.selected==TRUE){
        # if requested, correct for multiple testing
        # only the lower triangle of the matrix is corrected
        pvalue.matrix=correctPvalMatrix(pvalue.matrix)
      }
      #print("BH-corrected p-value matrix")
      #print(pvalue.matrix)
      # apply p-value threshold
      pvalue.matrix[pvalue.matrix>pval.T.selected]=NA
      # filter out edges not supported by the selected number of methods
      if(!is.na(method.num.T)){
        pvalue.matrix[score.matrix<method.num.T]=NA
      }
      sig.matrix=-1*log10(pvalue.matrix)
      sig.matrix[is.na(sig.matrix)]=0
      diag(sig.matrix)=0 # no self-loops

      # only the lower triangle will be used (this is important, because BH correction is only applied to the lower triangle of the p-value matrix)
      res.graph=graph_from_adjacency_matrix(sig.matrix,mode="lower",weighted=TRUE)
    }
  }

  adjacency.matrix=as_adj(res.graph, type="both", names=TRUE)
  #print(adjacency.matrix[2:10,])
  colors=c()
  nodecolors=c() # only needed for the bipartite case
  # assign signs as colors
  # graph is symmetric
  # igraph doc: The order of the vertices are preserved, i.e. the vertex corresponding to the first row will be vertex 0 in the graph, etc.
  for(i in 1:N){
    # color nodes of both types differently
    if(!is.null(metadata)){
      if(i>length(taxon.names)){
        nodecolors=c(nodecolors,"lightblue")
      }else{
        nodecolors=c(nodecolors,"salmon")
      }
    }
    for(j in 1:i){
      # skip diagonal
      if(i != j){
        # edge exists
        if(adjacency.matrix[i,j]!=0){
          if(sign.matrix[i,j]==1){
            colors=c(colors,"green")
          }else if(sign.matrix[i,j]==-1){
            colors=c(colors,"red")
          }else{
            colors=c(colors,"gray")
          }
        } # edge exists
      } # skip diagonal
    } # inner loop
  } # outer loop

  #print(length(colors))

  if(length(E(res.graph))==0){
    warning("The inferred network is empty.")
  }else{
    E(res.graph)$color=colors
    if(!is.null(metadata)){
      V(res.graph)$color=nodecolors
    }
    print(paste("Network has",length(E(res.graph)),"edges."))
    # do statistics
    num.pos=length(which(colors=="green"))
    num.neg=length(which(colors=="red"))
    num.uncertain=length(which(colors=="gray"))
    total=length(E(res.graph))
    one.percent=total/100
    pos.percent=num.pos/one.percent
    neg.percent=num.neg/one.percent
    uncertain.percent=num.uncertain/one.percent
    if(verbose == T){
      print(paste(num.pos,"positive edges",sep=" "))
      print(paste(num.neg,"negative edges",sep=" "))
      print(paste(num.uncertain,"conflicting edges",sep=" "))
      print(paste(pos.percent,"positive percent",sep=" "))
      print(paste(neg.percent,"negative percent",sep=" "))
      print(paste(uncertain.percent,"conflicting percent",sep=" "))
    }
  }
  # remove orphan nodes
  res.graph=delete.vertices(res.graph,igraph::degree(res.graph)==0) # degree is hidden by bnlearn
  return(res.graph)
}

# Compute row-wise associations and their p-values for the selected method. If thresholds are provided, apply them.
computeAssociations<-function(abundances, forbidden.combis=NULL, method="bray", permut=FALSE, bh=TRUE, T.down=NA, T.up=NA, pval.T=0.05, pval.cor=FALSE, renorm=FALSE, permutandboot=TRUE, iters=1000, verbose=FALSE, plot=FALSE, pseudocount=NA){
  print(paste("Renormalisation is",renorm))
  print(paste("P-values of correlations are computed with cor.test",pval.cor))
  print(paste("Permutations and bootstraps are both computed",permutandboot))

  N=nrow(abundances)
  #print(paste("p-value from cor.test?",pval.cor))
  correlations=c("spearman","pearson")
  #print(dim(forbidden.combis))

  sign.matrix = matrix(NA,nrow=N,ncol=N)

  ### compute p-values
  pvals = matrix(NA,nrow=N,ncol=N)
  if(permut == TRUE || permutandboot == TRUE || pval.cor == TRUE){
    for(i in 1:N){
      for(j in 1:i){
        # skip diagonal
        if(j!=i){
          # check whether combination is forbidden
          if(is.null(forbidden.combis) || !is.na(forbidden.combis[i,j])){
            # pval.cor takes precedence, but is only applied to correlations
            if(pval.cor == TRUE && method %in% correlations){
              pvals[i, j] = cor.test(abundances[i,],abundances[j,],method=method)$p.value
            }else{
              pvals[i, j] = getPval(abundances, i, j, method=method, N.rand=iters, renorm=renorm, permutandboot=permutandboot, verbose=verbose,  pseudocount=pseudocount)
              #print(pvals[i,j])
            }
          } # check forbidden combis
          pvals[j, i] = pvals[i, j]
        } # avoid diagonal
      } # inner loop
    } # outer loop
    # if requested, carry out Benjamini-Hochberg correction
    if(bh == TRUE){
      pvals=correctPvalMatrix(pvals)
    }
  }
  #print(pvals[!is.na(pvals)])

  # if provided, set forbidden combinations
  if(!is.null(forbidden.combis)){
    # scores have been computed previously for the thresholds or bipartite network and are re-used
    # sign matrix has been computed previously
    scores=forbidden.combis
  }else{
    scores=getScores(abundances,method=method,pseudocount=pseudocount)
    #print(scores[1,1:10])
    # compute sign matrix
    sign.matrix=computeSignMatrix(scores,method)
  }
  #print(scores)

  ### apply filter
  # if p-values were calculated, set all correlations with p-value above 0.05 to 0, set Bray Curtis scores to 0.5 and KLD scores to a value between the lower and upper threshold
  if(permut == TRUE || permutandboot == TRUE || pval.cor == TRUE){
    FILTER=pvals>pval.T
    pvals[FILTER]=NA # allows to set them to 0 later for adjacency conversion
  # apply thresholds on scores
  }else{
    if(!is.na(T.up) && !is.na(T.down)){
        print("up and down")
        # set all scores above lower threshold and below upper threshold to NA
        scores.temp=scores
        scores.temp[scores.temp>T.up]=Inf
        scores.temp[scores.temp<T.down]=Inf
        scores.temp[!is.infinite(scores.temp)]=NA
        scores[is.na(scores.temp)]=NA
    }else if(!is.na(T.up)){
        #print("up")
        scores[scores>T.up]=NA
    }else if(!is.na(T.down)){
        #print("down")
        scores[(scores>T.down)]=NA
    }
  }
  #print(scores)

  storage=scores
  # get lower triangle of score matrix
  scores=scores[lower.tri(scores)]
  # discard missing values for statistics and plotting
  scores=scores[!is.na(scores)]

  if(plot == T){
    values = scores
    what = "Score"
    if(permut == TRUE || permutandboot == TRUE || pval.cor==TRUE){
      values = pvals
      what = "P-value"
    }
    hist(values, main=paste(what," distribution for method ",method,sep=""), xlab=paste(what,"s",sep=""), ylab="Frequency")
  }

  # return score matrix (storage), pvals and sign matrix
  out=list(storage, pvals,sign.matrix)
  names(out)=c("scores","pvalues","signs")
  return(out)
}

# enforce bipartiteness given the separating row index and the matrix with forbidden combinations
# separator.index: first index where rows from metadata start
forbiddenCombisBipartite<-function(separator.index=NA,forbidden.combis=NULL){
  for(i in 1:nrow(forbidden.combis)){
    for(j in 1:nrow(forbidden.combis)){
      # edge is within node set 1
      if(i<separator.index && j<separator.index){
        forbidden.combis[i,j]=NA
      }
      # edge is within node set 2
      if(i>=separator.index && j>=separator.index){
        forbidden.combis[i,j]=NA
      }
    }
  }
  return(forbidden.combis)
}

# scores is a score matrix computed with the given method
computeSignMatrix<-function(scores, method="bray"){
  sign.matrix = matrix(NA,nrow=nrow(scores),ncol=ncol(scores))
  correlations=c("spearman","pearson")
  if(method %in% correlations){
    sign.matrix[scores>=0]=1 # copresence
    sign.matrix[scores<0]=-1 # mutual exclusion
  }else{
    if(method=="bray"){
      # Bray Curtis is bounded between 0 and 1, 1 being maximal dissimilarity
      sign.matrix[scores>0.5]=-1 # exclusion
      sign.matrix[scores<=0.5]=1 # copresence
    }else if(method=="kld"){
        mid.point=mean(range(scores))
        sign.matrix[scores>mid.point]=-1 # exclusion
        sign.matrix[scores<=mid.point]=1 # copresence
    }
  }
  return(sign.matrix)
}

# correct a pvalue matrix for multiple testing with Benjamini-Hochberg
# only the lower p-value matrix is corrected, but only the lower will be used
# to build the graph from the adjacency matrix
correctPvalMatrix<-function(pvals){
  pvec=pvals[lower.tri(pvals)]
  pvec=p.adjust(pvec,method="BH") # p.adjust keeps NA at the original place
  pvals[lower.tri(pvals)]=pvec
  return(pvals)
}

# get a symmetric row-wise association score matrix
getScores<-function(mat, method="spearman", pseudocount=NA){
  if(method == "spearman"){
    scores=cor(t(mat),method="spearman")
  }else if(method == "pearson"){
    scores=cor(t(mat),method="pearson")
  }else if(method == "bray"){
    scores = vegdist(mat, method="bray")
    scores=as.matrix(scores)
  }else if(method == "kld"){
    scores = computeKld(mat, pseudocount=pseudocount)
  }else{
    stop("Choose either spearman, pearson, kld or bray as a method.")
  }
  return(scores)
}


# Filter taxa in an abundance matrix
# Discard taxa with less than the given minimum number of occurrences.
# x taxon abundance matrix, rows are taxa, columns are samples
# minocc minimum occurrence (minimum number of samples with non-zero taxon abundance)
# keepSum If keepSum is true, the discarded rows are summed and the sum is added as a row with name: summed-nonfeat-rows
# return.filtered.indices if true, return an object with the filtered abundance matrix in mat and the indices of removed taxa in the original matrix in filtered.indices
# filtered abundance matrix
filterTaxonMatrix<-function(x, minocc=0, keepSum=FALSE, return.filtered.indices=FALSE){
  if(minocc==0){
    return(x)
  }else{
    toFilter=c()
    xcopy=x
    # convert into presence/absence matrix
    xcopy[xcopy>0]=1
    # sum for each taxon = number of occurrences across samples
    rowsums=apply(xcopy,1,sum)
    toFilter=which(rowsums<minocc)
    indices.tokeep=setdiff(c(1:nrow(x)),toFilter)
    #print(paste("Filtering",rownames(x)[toFilter]))
    if(keepSum==TRUE){
      filtered=x[toFilter,]
      x=x[indices.tokeep,]
      rownames=rownames(x)
      sums.filtered=apply(filtered,2,sum)
      x=rbind(x,sums.filtered)
      rownames=append(rownames,"summed-nonfeat-rows")
      rownames(x)=rownames
    }else{
      x=x[indices.tokeep,]
    }
    if(return.filtered.indices==TRUE){
      res=list(x,toFilter)
      names(res)=c("mat","filtered.indices")
      return(res)
    }else{
      return(x)
    }
  }
}


# Normalize a matrix
#
# Normalize a matrix column-wise by dividing each entry by its corresponding column sum.
#
# Columns summing to zero are removed by default.
#
# a matrix
# a normalized matrix
normalize<-function(x){
  # remove columns with only zeros from matrix, to avoid dividing by a zero
  colsums = apply(x,2,sum)
  for(i in 1:ncol(x)){
    x[,i]=x[,i]/colsums[i]
  }
  x
}


# Get p-value
# Get permuation-based p-value for association between two vectors.
#
# Compute the association between two vectors using the given method and
# compute its p-value using a permutation test. This method was adapted from R code by Fah Sahtirapongsasuti.

# matrix input matrix
# x.index index of first vector in the input matrix
# y.index index of second vector in the input matrix
# N.rand number of iterations used for the permutation test
# method similarity measure (supported measures are: "pearson", "spearman", "bray" and "kld")
# renorm renormalize after permutation
# permutandboot compute a bootstrap distribution in addition to the permutation distribution and
#                 return the p-value as the mean of the permutation distribution under the bootstrap distribution
# plot plot the histogram of the permutation and, if permutandboot is true, of both the permutation and the bootstrap distribution
# verbose print distribution properties and p-value
# pseudocount pseudocount used when computing KLD
#
# p-value of the association
#
# # Test
# data("ibd_taxa")
# data("ibd_lineages")
# ibd_genera=aggregateTaxa(ibd_taxa,lineages = ibd_lineages,taxon.level = "genus")
# getPval(matrix=ibd_genera,x.index=9,y.index=5,method="spearman",permutandboot=T, verbose=T,plot=T)
getPval = function(matrix, x.index, y.index, N.rand=1000, method="spearman", renorm=F, permutandboot=F, plot=F, verbose=F,  pseudocount=0.00000001) {
  x = as.numeric(matrix[x.index,])
  y = as.numeric(matrix[y.index,])
  lower.tail = FALSE
  #print(paste("Renorm applied:",renorm))
  if(method == "spearman"){
    this.sim = cor(x, y, use="complete.obs", method="spearman")
  }else if(method == "pearson"){
    this.sim = cor(x, y, use="complete.obs", method="pearson")
  }else if(method == "bray"){
    this.sim= vegdist(rbind(x,y),method="bray")
  }else if(method == "kld"){
    this.sim=get.kld(x,y, pseudocount = pseudocount)
  }else{
    stop("Select either spearman, pearson, kld or bray as method.")
  }
  rand.sim = rep(NA, N.rand)
  boot.sim = rep(NA, N.rand)
  for (i in 1:N.rand) {
    rand = sample(x, length(x))
    if(renorm == T){
      mat.copy=matrix
      mat.copy[x.index,]=rand
      mat.copy = normalize(mat.copy)
      rand = mat.copy[x.index,]
      y = mat.copy[y.index,]
    }
    if(method == "spearman"){
      rand.sim[i] = cor(rand, y, method="spearman", use="complete.obs")
    }else if(method == "pearson"){
      rand.sim[i] = cor(rand, y, method="pearson",use="complete.obs")
    }else if(method == "bray"){
      rand.sim[i] = vegdist(rbind(rand,y),method="bray")
    }else if(method == "kld"){
      rand.sim[i] = get.kld(rand,y, pseudocount = pseudocount)
    }
  }
  rand.sim = na.omit(rand.sim)
  if(plot == TRUE && permutandboot == FALSE){
    col1=rgb(0,0,1,1/3)
    col2=rgb(1,0,0,1/3)
    hist(rand.sim,col=col1)
    abline(v=mean(rand.sim),col="blue")
  }
  if(permutandboot){
    x=as.numeric(matrix[x.index,])
    y=as.numeric(matrix[y.index,])
    for (i in 1:N.rand) {
      rand.idx = sample(1:length(x),replace=TRUE)
      x.boot=x[rand.idx]
      y.boot=y[rand.idx]
      if(method == "spearman"){
        boot.sim[i] = cor(x.boot, y.boot, method="spearman", use="complete.obs")
      }else if(method == "pearson"){
        boot.sim[i] = cor(x.boot, y.boot, method="pearson",use="complete.obs")
      }else if(method == "bray"){
        boot.sim[i] = vegdist(rbind(x.boot,y.boot),method="bray")
      }else if(method == "kld"){
        boot.sim[i] = get.kld(x.boot,y.boot, pseudocount = pseudocount)
      }
     # print(boot.sim[i])
    }
    boot.sim = na.omit(boot.sim)
    if(plot == TRUE){
      compareDistribsPure(rand.sim,boot.sim, main="Permutation and bootstrap distributions", name1="Permutations", name2="Bootstraps", xlab=paste(method,"values",sep=" "))
    }
    # if we got enough non-NA permutation and bootstrap values, compute p-value
    if(length(rand.sim) > round(N.rand/3) && length(boot.sim) > round(N.rand/3)){
      # determine tail based on mean of permutation and bootstrap distributions
      if(mean(boot.sim)>mean(rand.sim)){
        lower.tail=TRUE
      }else{
        lower.tail=FALSE
      }
      pval = pnorm(mean(rand.sim),mean=mean(boot.sim),sd=sd(boot.sim), lower.tail=lower.tail)
    }else{
      pval = 0.5
    }
  }else{
    # if we got enough non-NA permutation values, compute p-value
    if(length(rand.sim) > round(N.rand/3)){
      # bray and kld are dissimilarities, so one-sided p-value needs to be computed from the lower tail (the smaller the better)
      if(method == "bray" || method == "kld"){
        lower.tail = TRUE
      }
      # look at left tail for negative correlations
      if(this.sim<0 && (method == "spearman" || method == "pearson")){
        lower.tail=TRUE
      }
      if (lower.tail) {
        pval = ((1+sum(this.sim > rand.sim))) / (1+length(rand.sim))
      } else {
        pval = ((1+sum(this.sim < rand.sim)) / (1+length(rand.sim)))
      }
    }else{
      pval = 0.5
    }
  }
  # set missing value (from constant vector) to intermediate non-significant p-value
  if(is.na(pval)){
    pval = 0.5
  }
  if(verbose == T){
    print(paste("p-value =",pval))
    print(paste("original score",this.sim))
    print(paste("mean of null distrib",mean(rand.sim)))
    print(paste("sd of null distrib",sd(rand.sim)))
    if(permutandboot == T){
      print(paste("mean of boot distrib",mean(boot.sim)))
      print(paste("sd of boot distrib",sd(boot.sim)))
    }
  }
  pval
}

# Compute KLD of a matrix row-wise
# Compute Kullback-Leibler dissimilarity
#
# Outputs a Kullback-Leibler dissimilarity (symmetrized divergence) matrix. Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#
# x matrix
# pseudocount pseudocount = this value is added to each zero
#
# kullback-leibler dissimilarity matrix
computeKld=function(x, pseudocount=0.00000001){
  # diagonal is zero
  kld=matrix(data=0,nrow=nrow(x),ncol=nrow(x))
  for(i in 1:nrow(x)){
    for(j in 1:i){
      kld[i,j]=get.kld(x[i,],x[j,], pseudocount=pseudocount)
      kld[j,i]=kld[i,j]
    }
  }
  kld
}

# Compute KLD between two vectors
# Compute Kullback-Leibler dissimilarity
#
# Outputs a Kullback-Leibler dissimilarity (symmetrized divergence) matrix. Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#
# x vector with non-negative numbers
# y vector with non-negative numbers
# pseudocount pseudocount = this value is added to each zero
#
# kullback-leibler dissimilarity value
get.kld=function(x,y, pseudocount=0.00000001){
  if(length(x) != length(y)){
    stop("The two vectors should have the same length!")
  }
  x[x==0]=pseudocount
  y[y==0]=pseudocount
  dis = 0
  x = x/sum(x)
  y = y/sum(y)
  for(i in 1:length(x)){
    if(!is.nan(x[i]) && !is.nan(y[i])){
      ratioxy = log(x[i]/y[i])
      ratioyx = log(y[i]/x[i])
      dis = x[i]*ratioxy+y[i]*ratioyx + dis
    }
  }
  dis
}
