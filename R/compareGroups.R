#' @title Compare diversity-related properties across groups
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param property richness (Chao1), evenness (Sheldon), alpha (alpha-diversity with Shannon index) or beta (beta-diversity)
#' @param method for beta-diversity, dissim (Bray Curtis dissimilarities) and DM (estimation of theta with the Dirichlet-Multnomial) are supported
#' @param groups group membership vector with as many entries as samples in abundances
#' @param plot.type pergroup (plot property per group; this is the only plot type available for richness, evenness, alpha and beta with DM) and intravsinter (plot property for all within-group samples versus all samples)
#' @param avg none, mean or median (average group property with the selected averaging method, does not work for method DM)
#' @param all include the beta diversity for all samples
#' @param noSameGroup for beta diversity with method dissim: only between-group pairs are used to calculate the beta diversity for all samples
#' @param rowNorm for beta diversity with method dissim: normalize the abundances row-wise
#' @param subsample for intravsinter plot type: subsample to have as many intra-group as inter-group pairs prior to performing a Wilcoxon test
#' @param xlab the x axis label
#' @param pvalViz if true and avg is set to none and plot.type is pergroup, significant Wilcoxon p-values are displayed on the box plot using function stat_compare_means in R package ggpubr
#'
#' @export
compareGroups<-function(abundances, property="beta", method="dissim", groups=c(), plot.type="pergroup", avg="none", all=FALSE, noSameGroup=TRUE, rowNorm=FALSE, subsample=TRUE, xlab="Groups", pvalViz=FALSE){
  if(plot.type=="intravsinter" && property!="beta"){
    stop("The intravsinter plot type is only supported for beta diversity!")
  }
  if(plot.type=="intravsinter" && method=="DM"){
    stop("The intravsinter plot type is not supported for method DM.")
  }
  if(length(groups)==0){
    warning("Empty group vector provided. All samples are assigned to the same group.")
    groups=rep("all",ncol(abundances))
  }
  if(property=="richness"){
    if(!is.whole.matrix(abundances)){
      vec=as.vector(abundances)
      scaling.factor=1/min(vec[vec!=0])
      warning(paste("Richness estimation requires counts; matrix is scaled to integers using factor:",scaling.factor))
      abundances=round(abundances*scaling.factor)
    }
  }
  if(property!="beta"){
    method=""
  }
  unique.groups=unique(groups)
  if(is.numeric(unique.groups)){
    unique.groups=sort(unique.groups)
  }
  # computing intra-group variability
  groups.with.theta=c()
  interdissim=c()
  intradissim=c()
  pergroupdissim=list()
  groupspecdissim=c()
  groupspecavgprops=c()
  main=""
  ylab=""

  if(rowNorm == TRUE){
    # normalize row-wise
    abundances=abundances/rowSums(abundances)
  }

  # loop groups
  for(group in unique.groups){
    groupspecdissim=c()
    print(paste("Processing group",group))
    group.indices=which(groups==group)
    group.data=abundances[,group.indices]
    if(property=="beta"){
      if(is.null(group.data)==FALSE && ncol(group.data) > 1){
        groups.with.theta=c(groups.with.theta,group)
        if(method == "dissim"){
          main="Beta diversity (Bray Curtis dissimilarities)"
          ylab="Bray Curtis dissimilarity"
          dissimMat=as.matrix(vegdist(t(group.data),method="bray"))
          dist=as.matrix(vegdist(t(group.data),method="bray"))
          pergroupdissim[[as.character(group)]]=dist[lower.tri(dist)] # paste("group",group,sep="")
          # average beta diversity
          if(avg!="none"){
            if(avg=="median"){
              groupspecavgprops=c(groupspecavgprops,median(groupspecdissim))
            }else if(avg=="mean"){
              groupspecavgprops=c(groupspecavgprops,mean(groupspecdissim))
            }
          }
        }else if(method == "DM"){
          main="Beta diversity (overdispersion)"
          ylab="Overdispersion"
          # dirmult: counts of alleles X (items=columns) vary across sub-populations (observations=rows), so transpose is necessary
          dm.fit=dirmult(t(group.data))
          intradissim=c(intradissim, dm.fit$theta)
        } # end beta diversity methods
      } # avoid patients with a single column
      else{
        warning(paste("Group",group,"has less than 2 samples: cannot compute beta diversity."))
      }
    }else if(property=="alpha" || property=="richness" || property=="evenness"){
      for(sample.index in 1:ncol(group.data)){
        if(property=="alpha"){
          main="Alpha diversity (Shannon index)"
          ylab="Shannon"
          # Shannon
          value=diversity(group.data[,sample.index])
        }else if(property=="richness"){
          main="Richness (Chao1 index)"
          ylab="Chao1"
          # Chao1
          value=estimateR(group.data[,sample.index])[2]
        }else if(property=="evenness"){
          main="Evenness (Sheldon index)"
          ylab="Sheldon"
          value=sheldon(group.data[,sample.index])
        }
        groupspecdissim=c(groupspecdissim,value)
        # average values over group members or keep all values
        if(avg=="mean"){
          groupspecavgprops=c(groupspecavgprops,mean(groupspecdissim))
        }else if(avg=="median"){
          groupspecavgprops=c(groupspecavgprops,median(groupspecdissim))
        }else if(avg=="none"){
          pergroupdissim[[as.character(group)]]=groupspecdissim
        }
      }
    }else{
      stop(paste("Property",property,"is not supported."))
    }
  } # end group loop

  # compute global property
  if(all){
    if(property=="beta"){
      if(method == "dissim"){
        dissimMat=as.matrix(vegdist(t(abundances),method="bray"))
        for(index1 in 1:(ncol(abundances)-1)){
          for(index2 in (index1+1):ncol(abundances)){
            group1=groups[index1]
            group2=groups[2]
            dissimVal = dissimMat[index1,index2]
            if(!noSameGroup || group1 != group2){
              interdissim=c(interdissim,dissimVal)
            }
          }
        }
      }else if(method == "DM"){
        dm.fit=dirmult(t(abundances))
        interdissim=dm.fit$theta
      }
    }else{
      warning("All is only supported for beta diversity.")
      all=FALSE
    }
  } # end all

  # plot
  if(plot.type=="pergroup"){
    if(method=="DM"){
      ntheta=length(intradissim)
      thetas=c(intradissim,interdissim)
      range=c(0,max(thetas))
      barplot(thetas,col=c(rep("green",ntheta),"red"),names.arg=c(groups.with.theta,"all"),xlab="Group",ylab="Theta",main="Estimated overdispersion",ylim=range)
    }else{
      if(avg!="none"){
        names(groupspecavgprops)=unique.groups
        if(all && property=="beta"){
          avgall=NA
          if(avg=="mean"){
            avgall=mean(interdissim)
          }else if(avg=="median"){
            avgall=median(interdissim)
          }
          groupspecavgprops=c(groupspecavgprops,avgall)
          names(groupspecavgprops)=c(unique.groups,"all")
        }
        range=c(0,max(groupspecavgprops))
        if(property=="beta"){
         range=c(0,1)
        }
        ylab=paste(avg,ylab)
        barplot(groupspecavgprops,xlab=xlab,ylab=ylab,main=main,ylim=range)
      }else{
        # richness, evenness and alpha-diversity that are not averaged as well as beta-diversity with method dissim are displayed as box plots
        if(all && property=="beta"){
          pergroupdissim[["all"]]=interdissim
        }
        mat=listToMat(pergroupdissim)
        ylim=c(0,max(mat,na.rm = TRUE))
        # display box plot with p-values
        if(pvalViz){
          df=as.data.frame(mat)
          # ggplot2 and reshape2 is imported by phyloseq, which is imported by seqgroup
          # stat_compare_means requires package ggpubr, also imported by seqgroup
          combinations=list()
          units=colnames(mat)
          for(index1 in 1:(length(units)-1)){
            for(index2 in (index1+1):length(units)){
              unit1=units[index1]
              unit2=units[index2]
              # two-sided, unpaired Wilcoxon test
              w.out=wilcox.test(mat[,index1],mat[,index2])
              if(w.out$p.value<0.05){
                #print(w.out$p.value)
                combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2)
              }
            }
          }
          print(paste("Number of significant differences in property",property,"across groups:",length(combinations)))
          variable="" # to avoid error message in package built
          df_melt <- melt(df)
          if(length(combinations)>0){
            # cannot set ylim, else p-values are not plotted correctly
            ggplot(df_melt, aes(variable,value))+geom_boxplot()+stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main)+geom_jitter(position = position_jitter(0.2))
          }else{
            ggplot(df_melt, aes(variable,value))+geom_boxplot()+stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main)+geom_jitter(position = position_jitter(0.2)) + ylim(ylim[1],ylim[2])
          }
        }else{
          boxplot(mat,ylab=ylab, main=main, xlab=xlab, notch=FALSE,ylim=ylim) # border=colors
          for(i in 1:ncol(mat)){
            points(rep(i,length(mat[,i])),mat[,i])
          }
        }
      }
    }
  }else if(plot.type=="intravsinter"){
    if(subsample){
      interdissimSub=interdissim[sample(1:length(intradissim))]
      mat=cbind(intradissim,interdissimSub)
    }else{
      dissimList=list(intradissim, interdissim)
      names(dissimList)=c("intra","inter")
      mat=listToMat(dissimList)
    }
    out=wilcox.test(mat[,1],mat[,2]) # unpaired, two-sided
    pval=round(out$p.value,4)
    colnames(mat)=c("intra-var","inter-var")
    min=min(mat)
    max=max(mat)
    maxD=max(max(hist(intradissim,plot=F)$counts),max(hist(interdissimSub,plot=F)$counts))
    cols=c(rgb(0,1,0,0.5),rgb(1,0,0,0.5))
    hist(intradissim,xlim=c(min,max),ylim=c(0,maxD), main=paste("P-value Wilcoxon",pval,sep=": "), xlab="Bray-Curtis dissimilarity", col=cols[1], breaks="FD")
    hist(interdissimSub,col=cols[2], breaks="FD", add=T)
    legend("topright", c("intra-var","inter-var"), lty = rep(1,length(cols)), col = cols, merge = TRUE, bg = "white", text.col="black")
  }else{
    stop(paste("Plot type",plot.type,"is not supported."))
  }
}

# Convert a list into a matrix. When a list entry
# has less values than the entry with the largest
# number of values, complete it with missing values.
listToMat<-function(groupprops=list()){
  # get group with the largest number of entries
  maxnum=0
  for(name in names(groupprops)){
    if(length(groupprops[[name]])>maxnum){
      maxnum=length(groupprops[[name]])
    }
  }
  mat=matrix(NA,nrow=maxnum,ncol=length(names(groupprops)))
  colnames(mat)=names(groupprops)
  counter=1
  for(name in names(groupprops)){
    for(index in 1:length(groupprops[[name]])){
      mat[index,counter]=groupprops[[name]][index]
    }
    counter=counter+1
  }
  return(mat)
}

# test whether matrix only contains integers.
is.whole.matrix<-function(x){
  if(length(which(apply(x,2,is.whole)==FALSE))>0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

# taken from: https://stat.ethz.ch/pipermail/r-help/2003-April/032471.html
is.whole <- function(a, tol = 1e-7) {
  is.eq <- function(x,y) {
    r <- all.equal(x,y, tol=tol)
    is.logical(r) && r
  }
  (is.numeric(a) && is.eq(a, floor(a))) ||
    (is.complex(a) && {ri <- c(Re(a),Im(a)); is.eq(ri, floor(ri))})
}

#' @title Compute evenness using Sheldon's index
#'
#' @description Sheldon's index is defined as \eqn{S=\frac{e^H}{N}}, where H is the Shannon diversity and N the species number.
#' It ranges from 0 to 1, where 1 signifies a perfectly even abundance distribution.
#'
#' @references A.L. Sheldon 1969. Equitability indices: dependence on the species count. Ecology, 50, 466-467.
#' @references C Heip 1974. A new index measuring evenness. J. mar. biol. Ass. UK 54, 555-557.
#'
#' @details Note that the N2N1 mode results in evenness smaller than 1 for equal taxon probabilities.
#'
#' @param x a vector of species abundances
#' @param correction whether or not to apply the correction described in Alatalo, Oikos 37, 199-204, 1981
#' @param N2N1 whether to compute Sheldon's evenness as the ratio of e raised to the power of H (H = Shannon diversity) and Simpson's diversity
#' @return Sheldon's evenness
##############################################################################
sheldon<-function(x, correction = TRUE, N2N1 = FALSE){
  H = diversity(x, index="shannon")
  if(N2N1){
    simpson = diversity(x, index="simpson")
    numerator = 1/simpson
    denominator = exp(1)^H
  }else{
    numerator = exp(1)^H
    denominator = specnumber(x)
  }
  if(correction){
    numerator = numerator - 1
    denominator = denominator - 1
  }
  S = numerator/denominator
  S
}
