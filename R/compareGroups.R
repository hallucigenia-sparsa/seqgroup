#' @title Compare diversity-related properties across groups
#' @description The function allows to compare richness, evenness, alpha diversity and beta diversity across groups, where
#' beta diversity can be assessed with both a Bray Curtis distribution or overdispersion from a Dirichlet-Multinomial distribution fit.
#' By default, groups are subsampled to the same sample number, for this reason the output of the same run can differ. To switch off this
#' behaviour, subsample can be set to false. Chao1 and Shannon diversity are computed with vegan.
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param property richness (Chao1), evenness (Sheldon), alpha (alpha-diversity with Shannon index) or beta (beta-diversity)
#' @param method for beta-diversity, dissim (Bray Curtis dissimilarities) and DM (estimation of theta with package dirmult) are supported
#' @param groups group membership vector with as many entries as samples in abundances
#' @param plot.type pergroup (plot property per group; this is the only plot type available for richness, evenness, alpha and beta with DM) and intravsinter (plot property for all within-group samples versus all samples)
#' @param avg none, mean or median (average group property with the selected averaging method, does not work for method DM)
#' @param all include the beta diversity for all samples
#' @param noSameGroup for beta diversity with method dissim: only between-group pairs are used to calculate the beta diversity for all samples
#' @param rowNorm for beta diversity with method dissim: normalize the abundances row-wise
#' @param subsample subsample groups randomly to have the same sample number in each group; for plot type intravsinter, it equalises in additon the number of pairs within and between groups
#' @param xlab the x axis label
#' @param pvalViz if true and avg is set to none and plot.type is pergroup, significant Wilcoxon p-values are displayed on the box plot using function stat_compare_means in R package ggpubr
#' @examples
#' data("ibd_taxa")
#' data("ibd_metadata")
#' groups=as.vector(ibd_metadata$Diagnosis)
#' compareGroups(ibd_taxa,groups=groups,property="alpha",pvalViz = TRUE)
#' @export
compareGroups<-function(abundances, property="beta", method="dissim", groups=c(), plot.type="pergroup", avg="none", all=FALSE, noSameGroup=TRUE, rowNorm=FALSE, subsample=TRUE, xlab="", pvalViz=FALSE){
  supported.properties=c("richness","evenness","alpha","beta")
  if(plot.type=="intravsinter" && property!="beta"){
    stop("The intravsinter plot type is only supported for beta diversity!")
  }
  if(plot.type=="intravsinter" && method=="DM"){
    stop("The intravsinter plot type is not supported for method DM.")
  }
  if(length(groups)==0){
    if(plot.type=="intravsinter"){
      stop("For plot type intravsinter, a vector of group memberships is required.")
    }else{
      warning("Empty group vector provided. All samples are assigned to the same group.")
      groups=rep("all",ncol(abundances))
    }
  }
  if(plot.type=="intravsinter"){
      all=TRUE
      property="beta"
  }
  if(method %in% supported.properties){
    stop(paste("The method parameter specifies the way in which beta diversity is computed. Please use the property parameter."))
  }
  scaleMatrix=FALSE
  if(property=="richness"){
    if(!is.Count.Matrix(abundances)){
      warning("Richness estimation requires counts. Data were scaled to counts.")
      scaleMatrix=TRUE
    }
  }
  if(property=="beta" && method=="DM"){
    if(!is.Count.Matrix(abundances)){
      warning("The Dirichlet-Multinomial distribution can only be estimated from count data. Data were scaled to counts.")
      scaleMatrix=TRUE
    }
  }
  if(scaleMatrix){
    vec=as.vector(abundances)
    scaling.factor=1/min(vec[vec!=0])
    if(scaling.factor>10000){
      scaling.factor=10000
    }
    print(paste("Abundances have been scaled with factor:",scaling.factor))
    abundances=round(abundances*scaling.factor)
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

  if(xlab==""){
    if(plot.type=="intravsinter"){
      xlab="Bray-Curtis dissimilarity"
    }else{
      xlab="Groups"
    }
  }

  if(rowNorm == TRUE){
    # normalize row-wise
    abundances=abundances/rowSums(abundances)
  }

  constrainSampleNum=FALSE
  sampleNum=NA

  if(subsample==TRUE){
    constrainSampleNum=TRUE
    sampleNum=min(table(groups))
    print(paste("Constraining sample number randomly to the same minimal group sample number of",sampleNum))
  }

  # loop groups
  for(group in unique.groups){
    groupspecdissim=c()
    print(paste("Processing group",group))
    group.indices=which(groups==group)
    if(constrainSampleNum){
      group.indices=sample(group.indices)[1:sampleNum]
    }
    group.data=abundances[,group.indices]
    if(property=="beta"){
      if(is.null(group.data)==FALSE && ncol(group.data) > 1){
        groups.with.theta=c(groups.with.theta,group)
        if(method == "dissim"){
          main="Beta diversity (Bray Curtis dissimilarities)"
          ylab="Bray Curtis dissimilarity"
          dissimMat=as.matrix(vegdist(t(group.data),method="bray"))
          dist=as.matrix(vegdist(t(group.data),method="bray"))
          pergroupdissim[[as.character(group)]]=dist[lower.tri(dist)]
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
          value=vegan::diversity(group.data[,sample.index])
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
      thetas=intradissim
      print(interdissim)
      names.arg=groups.with.theta
      if(all==TRUE){
        names.arg=c(groups.with.theta,"all")
        thetas=c(intradissim,interdissim)
      }
      range=c(0,max(thetas))
      barplot(thetas,names.arg=names.arg,xlab=xlab,ylab="Theta",main="Estimated overdispersion",ylim=range)
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
            ggplot(df_melt, aes(variable,value))+geom_boxplot()+ggpubr::stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main)+geom_jitter(position = position_jitter(0.2))
          }else{
            ggplot(df_melt, aes(variable,value))+geom_boxplot()+ggpubr::stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main)+geom_jitter(position = position_jitter(0.2)) + ylim(ylim[1],ylim[2])
          }
        }else{
          #print(dim(mat))
          #print(colnames(mat))
          boxplot(mat,ylab=ylab, main=main, xlab=xlab, notch=FALSE,ylim=ylim) # border=colors
          for(i in 1:ncol(mat)){
            points(rep(i,length(mat[,i])),mat[,i])
          }
        }
      }
    }
  }else if(plot.type=="intravsinter"){
    intradissim=as.vector(listToMat(pergroupdissim))
    print(paste("Number of inter-group values:",length(interdissim)))
    print(paste("Number of intra-group values:",length(intradissim)))
    if(subsample){
      print("Number of inter and intra-group values is equalised.")
      if(length(intradissim)<length(interdissim)){
        mat=cbind(intradissim,interdissim[sample(1:length(intradissim))])
      }else{
        mat=cbind(intradissim[sample(1:length(interdissim))],interdissim)
      }
    }else{
      dissimList=list(intradissim, interdissim)
      names(dissimList)=c("intra","inter") # names are required here
      mat=listToMat(dissimList)
    }
    colnames(mat)=c("intra","inter")
    out=wilcox.test(mat[,1],mat[,2],na.rm=TRUE) # unpaired, two-sided
    pval=round(out$p.value,4)
    min=min(mat,na.rm =TRUE)
    max=max(mat,na.rm =TRUE)
    maxD=max(max(hist(mat[,1],plot=F)$counts, na.rm=TRUE),max(hist(mat[,2],plot=F)$counts, na.rm=TRUE))
    #print(min)
    #print(max)
    #print(maxD)
    cols=c(rgb(0,1,0,0.5),rgb(1,0,0,0.5))
    hist(mat[,1],xlim=c(min,max),ylim=c(0,maxD), main=paste("P-value Wilcoxon",pval,sep=": "), xlab=xlab, col=cols[1], breaks="FD")
    hist(mat[,2],col=cols[2], breaks="FD", add=T)
    legend("topleft", colnames(mat), lty = rep(1,length(cols)), col = cols, merge = TRUE, bg = "white", text.col="black")
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
  #print(names(groupprops))
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

# Compute evenness using Sheldon's index
#
# Sheldon's index is defined as \eqn{S=\frac{e^H}{N}}, where H is the Shannon diversity and N the species number.
# It ranges from 0 to 1, where 1 signifies a perfectly even abundance distribution.
#
#  A.L. Sheldon 1969. Equitability indices: dependence on the species count. Ecology, 50, 466-467.
# C Heip 1974. A new index measuring evenness. J. mar. biol. Ass. UK 54, 555-557.
#
#  Note that the N2N1 mode results in evenness smaller than 1 for equal taxon probabilities.
#
# x a vector of species abundances
# correction whether or not to apply the correction described in Alatalo, Oikos 37, 199-204, 1981
# N2N1 whether to compute Sheldon's evenness as the ratio of e raised to the power of H (H = Shannon diversity) and Simpson's diversity
# Sheldon's evenness
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
