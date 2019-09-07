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
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param methods network construction methods, values can be combinations of: "pearson", "spearman", "kld" or "bray"
#' @param T.up upper threshold for scores (when more than one network construction method is provided, T.up is ignored)
#' @param T.down lower threshold for scores (when more than one network construction method is provided, T.down is ignored)
#' @param method.num.T threshold on method number (only used when more than one method is provided)
#' @param pval.T threshold on p-value (only used when permut is true); if several methods are provided, only applied after merge
#' @param init.edge.num the number of top and bottom initial edges (init.edge.num overrides T.up/T.down, set to NA to respect T.up/T.down for a single method)
#' @param min.occ only keep rows with at least the given number of non-zero values (carried out before network construction)
#' @param keep.filtered sum all filtered rows and add the sum vector as additional row
#' @param norm normalize matrix (carrried out after filtering)
#' @param stand.rows standardize rows by dividing each entry by its corresponding row sum, applied after normalization
#' @param pval.cor compute p-values of correlations with cor.test (only valid for correlations; takes precedence over permut and/or permutandboot with or without renorm)
#' @param permut compute p-values on edges with a permutation test
#' @param renorm compute p-values with a permutation test, using renormalization (only applied to correlations)
#' @param permutandboot compute p-values from both permutation (with or without renorm) and bootstrap distribution
#' @param iters number of iterations for the permutation test
#' @param bh multiple-test-correct using Benjamini-Hochberg; if several methods are provided bh is applied to merged p-value
#' @param pseudocount count added to zeros prior to taking logarithm (for KLD, p-value merge and significance)
#' @param plot plot score or, if permut or pval.cor is true, p-value distribution
#' @param verbose print the number of positive and negative edges and, if permut is true, details of p-value computation
#' @return igraph object with absolute association strengths, number of supporting methods or, if permut or pval.cor is true, significances (-1*log10(pval)) as edge weights
#' @examples
#' data("ibd_taxa")
#' data("ibd_lineages")
#' ibd_genera=aggregateTaxa(ibd_taxa,lineages = ibd_lineages,taxon.level = "genus")
#' min.occ=nrow(ibd_genera)/3
#' # p-values for the 50 strongest positive and 50 strongest negative Spearman correlations
#' plot(barebonesCoNet(ibd_genera,methods="spearman",init.edge.num=50,min.occ=min.occ,pval.cor=TRUE))
#' # combine Bray Curtis and Spearman without computing p-values
#' plot(barebonesCoNet(ibd_genera,methods=c("spearman","bray"),init.edge.num=50,min.occ = min.occ))
#' @export
barebonesCoNet<-function(abundances, methods=c("spearman","kld"), T.up=NA, T.down=NA, method.num.T=2, pval.T=0.05, init.edge.num=max(2,round(sqrt(nrow(abundances)))), min.occ=0, keep.filtered=TRUE, norm=FALSE, stand.rows=FALSE, pval.cor=FALSE, permut=FALSE, renorm=FALSE, permutandboot=FALSE, iters=100, bh=TRUE, pseudocount=0.00000000001, plot=FALSE, verbose=FALSE){

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

  if(permut==TRUE && iters < 1){
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

  ### Preprocessing
  abundances = filterTaxonMatrix(abundances,keepSum = keep.filtered, minocc=min.occ)
  N=nrow(abundances)
  # normalize matrix
  if(norm == TRUE){
    abundances=normalize(abundances)
  }
  # normalize row-wise
  if(stand.rows == TRUE){
    abundances = t(normalize(t(abundances)))
  }

  ### Network construction
  # loop selected methods
  for(method in methods){
    sign.matrix=matrix(NA,nrow=N,ncol=N)
    renorm.selected=renorm
    print(paste("Processing method",method))
    if(!is.na(init.edge.num) && init.edge.num>0){
      #  get score matrix
      scores=getScores(abundances,method=method, pseudocount=pseudocount)
      #print(length(as.vector(scores[lower.tri(scores)])))
      # sort score vector representing lower triangle of score matrix in ascending order
      scorevec=sort(as.vector(scores[lower.tri(scores)]),decreasing = FALSE)
      #print(scorevec[1:10])
      # determine lower threshold
      T.down=scorevec[init.edge.num]
      # determine upper threshold
      T.up=scorevec[(length(scorevec)-init.edge.num)]
      print(paste("Lower threshold for initial edge number (",init.edge.num,"): ",T.down,sep=""))
      print(paste("Upper threshold for initial edge number (",init.edge.num,"): ",T.up,sep=""))
      # compute sign matrix given thresholds
      if(method %in% correlations){
        sign.matrix[scores>=T.up]=1 # copresence
        sign.matrix[scores<=T.down]=-1 # mutual exclusion
      }else{
        # dissimilarity: scores above T.up represent exclusion
        sign.matrix[scores>=T.up]=-1 # exclusion
        sign.matrix[scores<=T.down]=1 # copresence
      }
      sign.matrix.list[[method]]=sign.matrix
      # set scores to avoid to NA
      forbidden.combis=scores
      forbidden.combis[forbidden.combis>T.up]=Inf
      forbidden.combis[forbidden.combis<T.down]=Inf
      forbidden.combis[!is.infinite(forbidden.combis)]=NA
      #forbidden.indices=which(forbidden.combis>T.down && forbidden.combis<T.up,arr.ind = TRUE)
      scores[is.na(forbidden.combis)]=NA
      forbidden.combis=scores
      #print(scores)
      # thresholds are not applied when initial edge number is indicated
      T.down=NA
      T.up=NA
      if(renorm == TRUE & (method == "kld" || method == "bray")){
        renorm.selected=FALSE # kld and bray are compositionality-robust
      }
    }
    res=computeAssociations(abundances,method=method, forbidden.combis = forbidden.combis, pval.T = pval.T, bh=bh, T.down=T.down, T.up=T.up, pval.cor = pval.cor, renorm=renorm.selected, permut=permut, permutandboot = permutandboot, verbose=verbose, plot=plot, iters=iters, pseudocount=pseudocount)
    resultList[[method]]=res
  }

  print(paste("Associations computed for",N,"taxa"))

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
    if(permut==TRUE || pval.cor==TRUE){
      # for a single method, multiple-testing correction was carried out previously
      # avoid taking the logarithm of zero p-value, else we loose them
      res$pvalues[res$pvalues==0]=pseudocount
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
       res$pvalues[res$pvalues==0]=pseudocount
       logpvalue.matrix=log(res$pvalues)
       # after log, set missing values to 0, so they do not set the entire calculation to NA (method is not counted)
       logpvalue.matrix[is.na(logpvalue.matrix)]=0
       # multiply logarithm of individual p-values
       if(permut==TRUE || pval.cor==TRUE){
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
    if(!permut && !pval.cor){
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
  # assign signs as colors
  # graph is symmetric
  # igraph doc: The order of the vertices are preserved, i.e. the vertex corresponding to the first row will be vertex 0 in the graph, etc.
  for(i in 1:N){
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
  res.graph=delete.vertices(res.graph,degree(res.graph)==0)
  return(res.graph)
}

# Compute row-wise associations and their p-values using the selected methods. If thresholds are provided, apply them.
computeAssociations<-function(abundances, forbidden.combis=NULL, method="bray", permut=FALSE, bh=TRUE, T.down=NA, T.up=NA, pval.T=0.05, pval.cor=FALSE, renorm=FALSE, permutandboot=TRUE, iters=1000, verbose=FALSE, plot=FALSE, pseudocount=NA){
  N=nrow(abundances)
  #print(paste("p-value from cor.test?",pval.cor))
  correlations=c("spearman","pearson")
  print(dim(forbidden.combis))

  ### compute p-values
  pvals = matrix(NA,nrow=N,ncol=N)
  if(permut == TRUE || pval.cor == TRUE){
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
    # scores have been computed previously for the thresholds and are re-used
    scores=forbidden.combis
  }else{
    scores=getScores(abundances,method=method,pseudocount=pseudocount)
    #print(scores[1,1:10])
  }
  #print(scores)

  ### apply filter
  # if p-values were calculated, set all correlations with p-value above 0.05 to 0, set Bray Curtis scores to 0.5 and KLD scores to a value between the lower and upper threshold
  if(permut == TRUE || pval.cor == TRUE){
    FILTER=pvals>pval.T
    pvals[FILTER]=NA # allows to set them to 0 later for adjacency conversion
  # apply thresholds on scores
  }else{
    if(!is.na(T.up) && !is.na(T.down)){
        #print("up and down")
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
    if(permut == TRUE || pval.cor==TRUE){
      values = pvals
      what = "P-value"
    }
    hist(values, main=paste(what," distribution for method ",method,sep=""), xlab=paste(what,"s",sep=""), ylab="Frequency")
  }

  # return score matrix (storage) and pvals matrix
  out=list(storage, pvals)
  names(out)=c("scores","pvalues")
  return(out)
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


#' @title Filter taxa in an abundance matrix
#' @description Discard taxa with less than the given minimum number of occurrences.
#' @param x taxon abundance matrix, rows are taxa, columns are samples
#' @param minocc minimum occurrence (minimum number of samples with non-zero taxon abundance)
#' @param keepSum If keepSum is true, the discarded rows are summed and the sum is added as a row with name: summed-nonfeat-rows
#' @param return.filtered.indices if true, return an object with the filtered abundance matrix in mat and the indices of removed taxa in the original matrix in filtered.indices
#' @return filtered abundance matrix
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


#' @title Normalize a matrix
#'
#' @description Normalize a matrix column-wise by dividing each entry by its corresponding column sum.
#'
#' @details Columns summing to zero are removed by default.
#'
#' @param x a matrix
#' @return a normalized matrix
normalize<-function(x){
  # remove columns with only zeros from matrix, to avoid dividing by a zero
  colsums = apply(x,2,sum)
  for(i in 1:ncol(x)){
    x[,i]=x[,i]/colsums[i]
  }
  x
}


#' @title Get p-value
#' @description Get permuation-based p-value for association between two vectors.
#'
#' @details # Compute the association between two vectors using the given method and
#' compute its p-value using a permutation test. This method was adapted from R code by Fah Sahtirapongsasuti.
#' This method was adapted from CCREPE: http://huttenhower.sph.harvard.edu/ccrepe.
#' Emma Schwager et al Detecting statistically significant associtations between sparse and high dimensional compositional data. (In progress).
#'
#' @param matrix input matrix
#' @param x.index index of first vector in the input matrix
#' @param y.index index of second vector in the input matrix
#' @param N.rand number of iterations used for the permutation test
#' @param method similarity measure (supported measures are: "pearson", "spearman", "bray" and "kld")
#' @param renorm renormalize after permutation
#' @param permutandboot compute a bootstrap distribution in addition to the permutation distribution and
#                 return the p-value as the mean of the permutation distribution under the bootstrap distribution
#' @param plot plot the histogram of the permutation and, if permutandboot is true, of the bootstrap distribution
#' @param verbose print distribution properties and p-value
#' @param pseudocount pseudocount used when computing KLD
#'
#' @return p-value of the association
getPval = function(matrix, x.index, y.index, N.rand=1000, method="spearman", renorm=F, permutandboot=F, plot=F, verbose=F,  pseudocount=NA) {
  x = matrix[x.index,]
  y = matrix[y.index,]
  lower.tail = TRUE
  #print(paste("Renorm applied:",renorm))
  # bray and kld are dissimilarities, so one-sided p-value needs to be computed from the upper tail
  if(method == "bray" || method == "kld"){
    lower.tail = FALSE
  }
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
  if(plot == T){
    col1=rgb(0,0,1,1/3)
    col2=rgb(1,0,0,1/3)
    hist(rand.sim,col=col1)
    abline(v=mean(rand.sim),col="blue")
  }
  if(permutandboot){
    x=matrix[x.index,]
    y=matrix[y.index,]
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
    }
    boot.sim = na.omit(boot.sim)
    if(plot == T){
      hist(boot.sim,col=col2,add=T)
      abline(v=mean(boot.sim),col="red")
      legend(x="topleft", c("Permut","Boot"), bg="white",col=c(col1,col2),lty=rep(1,2),merge=T)
    }
    # if we got enough non-NA permutation and bootstrap values, compute p-value
    if(length(rand.sim) > round(N.rand/3) && length(boot.sim) > round(N.rand/3)){
      pval = pnorm(mean(rand.sim),mean=mean(boot.sim),sd=sd(boot.sim), lower.tail=lower.tail)
    }else{
      pval = 0.5
    }
  }else{
    # if we got enough non-NA permutation values, compute p-value
    if(length(rand.sim) > round(N.rand/3)){
      if (lower.tail) {
        pval = (sum(this.sim > rand.sim) / length(rand.sim))
      } else {
        pval = (sum(this.sim < rand.sim) / length(rand.sim))
      }
    }else{
      pval = 0.5
    }
  }
  # set missing value (from constant vector) to intermediate p-value (worst possible p-value in this context)
  if(is.na(pval)){
    pval = 0.5
  }
  # p-values are one-sided, so high p-values signal mutual exclusion and are converted into low ones
  if(pval > 0.5){
    pval = 1 - pval
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

#' @title Compute KLD of a matrix row-wise
#' @description Compute Kullback-Leibler dissimilarity
#'
#' @details Outputs a Kullback-Leibler dissimilarity (symmetrized divergence) matrix. Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#'
#' @param x matrix
#' @param pseudocount pseudocount = this value is added to each zero
#
#' @return kullback-leibler dissimilarity matrix
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

#' @title Compute KLD between two vectors
#' @description Compute Kullback-Leibler dissimilarity
#'
#' @details Outputs a Kullback-Leibler dissimilarity (symmetrized divergence) matrix. Note:
# Equation: D(x,y) = SUM(x_i*log(x_i/y_i) + y_i*log(y_i/x_i))
# taken from "Caution! Compositions! Can constraints on omics data lead analyses astray?"
# David Lovell et al., Report Number EP10994
#'
#' @param x vector with non-negative numbers
#' @param y vector with non-negative numbers
#' @param pseudocount pseudocount = this value is added to each zero
#
#' @return kullback-leibler dissimilarity value
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
