#' @title Compare diversity-related properties across groups
#' @description The function allows to compare richness, evenness, alpha diversity and beta diversity across groups, where
#' beta diversity can be assessed with both a Bray Curtis distribution or overdispersion from a Dirichlet-Multinomial distribution fit.
#' By default, groups are subsampled to the same sample number, for this reason the output of the same run can differ. To switch off this
#' behaviour, subsample can be set to false. Chao1 and Shannon diversity are computed with vegan.
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param property richness (Chao1), evenness (Sheldon), alpha (alpha-diversity with Shannon index) or beta (beta-diversity)
#' @param method for beta-diversity, dissim (Bray Curtis dissimilarities) and DM (estimation of theta with package dirmult) are supported
#' @param groups group membership vector with as many entries as samples in abundances
#' @param colors optional color vector with as many entries as (unique) groups (only for plot.type pergroup)
#' @param colorMap optional color map with a color for a metadata value, required for color legend (only for plot.type pergroup)
#' @param plot.type pergroup (plot property per group; this is the only plot type available for richness, evenness, alpha and beta with DM) and intravsinter (plot property for all within-group samples versus all samples)
#' @param avg none, mean or median (average group property with the selected averaging method, does not work for method DM)
#' @param all include the beta diversity for all samples
#' @param noSameGroup for beta diversity with method dissim: only between-group pairs are used to calculate the beta diversity for all samples
#' @param rowNorm for beta diversity with method dissim: normalize the abundances row-wise
#' @param subsample subsample groups randomly to have the same sample number in each group; for plot type intravsinter, it equalises in additon the number of pairs within and between groups
#' @param xlab the x axis label
#' @param pvalViz if true and avg is set to none and plot.type is pergroup, significant Wilcoxon p-values are displayed on the violin plot using function stat_compare_means in R package ggpubr
#' @param returnValues if true, return properties
#' @return if returnValues is true, properties are returned
#' @examples
#' data("ibd_taxa")
#' data("ibd_metadata")
#' groups=as.vector(ibd_metadata$Diagnosis)
#' compareGroups(ibd_taxa,groups=groups,property="alpha",pvalViz = TRUE)
#' data("floresgut_taxa")
#' data("floresgut_metadata")
#' groups=floresgut_metadata$host_subject_id
#' diet.col=makeColorsGivenMetadata(metadataItem=floresgut_metadata$diet,groups=groups)
#' compareGroups(floresgut_taxa,groups=groups,colors=diet.col$colors,colorMap=diet.col$colormap)
#' age=floresgut_metadata$age
#' age.col=makeColorsGivenMetadata(metadataItem=age,groups=groups,binning="equaldist",numBins=10)
#' col=age.col$colors
#' compareGroups(floresgut_taxa,groups=groups,property="alpha",colors=col,colorMap=age.col$colormap)
#' @export
compareGroups<-function(abundances, property="beta", method="dissim", groups=c(), colors=c(), colorMap=NULL, plot.type="pergroup", avg="none", all=FALSE, noSameGroup=TRUE, rowNorm=FALSE, subsample=TRUE, xlab="", pvalViz=FALSE, returnValues=FALSE){
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

  las=1
  if(length(unique(groups))>10){
    las=2 # vertical x axis labels
  }

  if(length(colors)>0 && length(groups)>0){
    if(length(colors)!=length(unique(groups))){
      #print(length(unique(groups)))
      #print(length(colors))
      stop("Please provide as many colors as there are groups!")
    }
  }

  if(length(colors)==0 && length(groups)>0){
    colors=rep("white",length(groups))
  }else if(length(colors)==0 && length(groups)==0){
    # default color if no groups are provided
    colors=c("black")
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

  # return
  res=NULL

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
    #print(paste("Processing group",group))
    group.indices=which(groups==group)
    if(constrainSampleNum){
      group.indices=sample(group.indices)[1:sampleNum]
    }
    group.data=abundances[,group.indices]
    if(length(group.indices)==1){
      group.data=as.matrix(group.data)
    }
    if(property=="beta"){
      if(is.null(group.data)==FALSE && ncol(group.data) > 0){
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

  if(las==2){
    xlab=""
  }
  #if(!is.null(colorMap)){
    # TODO to make space in par for the legend and shift it there, no longer needed for violin plots, but still required for barplots
  #}

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
      res=thetas
      compareGroupsPlot(mat=thetas, labels=names.arg, xlab=xlab, ylab="Theta", main="Estimated overdispersion", range=range, plotType = "bar",colors=colors, las=las, colorMap=colorMap)
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
        res=groupspecavgprops
        compareGroupsPlot(mat=groupspecavgprops, xlab=xlab, ylab=ylab, main=main, range=range, plotType = "bar",colors=colors, las=las, colorMap=colorMap)
      }else{
        # richness, evenness and alpha-diversity that are not averaged as well as beta-diversity with method dissim are displayed as box plots
        if(all && property=="beta"){
          pergroupdissim[["all"]]=interdissim
        }
        mat=listToMat(pergroupdissim)
        ylim=c(0,max(mat,na.rm = TRUE))
        res=mat
        # display box plot with p-values
        if(pvalViz){
          #df=as.data.frame(mat)
          combinations=list()
          units=colnames(mat)
          for(index1 in 1:(length(units)-1)){
            for(index2 in (index1+1):length(units)){
              unit1=units[index1]
              unit2=units[index2]
              # two-sided, unpaired Wilcoxon test
              val1=unique(mat[,index1])
              val2=unique(mat[,index2])
              isProblem=FALSE
              if((length(val1)==1 && is.na(val1)) || (length(val2)==1 && is.na(val2))){
                isProblem=TRUE
              }
              # avoid result only consisting of missing values (e.g. if group has only 1 sample)
              if(!isProblem){
                w.out=wilcox.test(mat[,index1],mat[,index2])
              }else{
                w.out$p.value=0.5
              }
              if(w.out$p.value<0.05){
                #print(w.out$p.value)
                combinations[[paste(unit1,unit2,sep="")]]=c(unit1,unit2)
              }
            }
          }
          print(paste("Number of significant differences in property",property,"across groups:",length(combinations)))
          compareGroupsPlot(mat=mat, groups=groups, colors=colors, combinations=combinations, xlab=xlab, ylab=ylab, main=main, plotType = "pvalues")
        }else{
          # alternative plot type: box, violin
          compareGroupsPlot(mat=mat, xlab=xlab, ylab=ylab, groups=groups, main=main, colors=colors, las=las, colorMap=colorMap, plotType = "violin")
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
    res=mat
    hist(mat[,1],xlim=c(min,max),ylim=c(0,maxD), main=paste("P-value Wilcoxon",pval,sep=": "), xlab=xlab, col=cols[1], breaks="FD")
    hist(mat[,2],col=cols[2], breaks="FD", add=T)
    legend("topleft", colnames(mat), lty = rep(1,length(cols)), col = cols, merge = TRUE, bg = "white", text.col="black")
  }else{
    stop(paste("Plot type",plot.type,"is not supported."))
  }

  if(returnValues){
    return(res)
  }
}

# Internal function for group-wise property plotting.
# mat: values to be plotted
# xlab: label for x axis
# ylab: label for y axis
# range: y axis range
# main: title of the plot
# groups: group membership vector
# labels: bar plot labels
# colors: bar/box color
# combinations: for pvalue plots
# las: orientation of axis labels
# plotType: box, bar, pvalues or violin
# colorMap: color map for a metadata item, adds color legend
compareGroupsPlot<-function(mat, xlab="", ylab="", range=c(), main="", groups=c(), labels=c(), colors=c(), combinations=c(), las=1, plotType="box", colorMap=NULL){

  legend.colors=c()
  if(!is.null(colorMap)){
    for(entry in names(colorMap)){
      legend.colors=c(legend.colors,colorMap[[entry]])
    }
  }

  if(plotType=="box"){
    # varwidth does not work, notch generates a warning message
    boxplot(mat,ylab=ylab, main=main, xlab=xlab, notch=FALSE,ylim=range, col = colors, cex.names = 0.9, las=las) # border=colors
    for(i in 1:ncol(mat)){
        points(rep(i,length(mat[,i])),mat[,i])
    }
  # make a plot with ggplot2
  }else if(plotType=="violin" || plotType=="pvalues"){
    # ggplot2 and reshape2 are imported by phyloseq, which is imported by seqgroup
    # stat_compare_means requires package ggpubr, also imported by seqgroup
    # to avoid error messages in package built
    variable=""
    value=""
    color=""
    df=as.data.frame(mat)
    df_melt <- melt(df)
    # if a color map is given, add color as third factor
    if(!is.null(colorMap)){
      # melt color map into two vectors
      legend.names=names(colorMap)
      legend.colors=c()
      for(name in names(colorMap)){
        legend.colors=c(legend.colors,colorMap[[name]])
      }
      #print(colors)
      #print(legend.names)
      #print(legend.colors)
      df_melt_extended=matrix(nrow=nrow(df_melt),ncol=3)
      unique.groups=unique(groups)
      #print(unique.groups)
      # loop rows of df_melt object and add color
      for(i in 1:nrow(df_melt)){
        group=as.character(df_melt[i,1])
        group.index=which(unique.groups==group)
        df_melt_extended[i,1]=group
        df_melt_extended[i,2]=as.numeric(as.character(df_melt[i,2]))
        # find the name of the metadata item in the color map
        color.name.index=which(legend.colors==as.character(colors[group.index]))
        df_melt_extended[i,3]=as.character(legend.names[color.name.index])
      }
      colnames(df_melt_extended)=c("variable","value","color")
      df_melt_extended=as.data.frame(df_melt_extended)
      # make sure value is not treated as a factor
      df_melt_extended$value=as.numeric(as.character(df_melt_extended$value))
    }
    if(plotType=="pvalues"){
      # cannot set ylim, else p-values are not plotted correctly
      if(!is.null(colorMap)){
        p <- ggplot(df_melt_extended, aes(variable, value, fill=color))+geom_violin()+ ggpubr::stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # 90 degree x axis test
        plot(p)
      }else{
        # geom_jitter(position = position_jitter(0.2))
        p <- ggplot(df_melt, aes(variable, value, fill=variable))+ guides(fill=FALSE) +geom_violin(show.legend=FALSE)+ scale_fill_manual(values=colors) + ggpubr::stat_compare_means(comparisons=combinations)+xlab(xlab)+ylab(ylab)+ggtitle(main) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # 90 degree x axis test
        plot(p)
      }
    }else{
      if(!is.null(colorMap)){
        p <- ggplot(df_melt_extended, aes(variable,value, fill=color)) + geom_violin()  +  xlab(xlab) + ylab(ylab) + ggtitle(main) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }else{
        # to remove legend: http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
        p <- ggplot(df_melt, aes(variable,value, fill=variable)) + guides(fill=FALSE) + geom_violin()  + scale_fill_manual(values=colors) + xlab(xlab) + ylab(ylab) + ggtitle(main) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        plot(p)
      }
    }
  }else if(plotType=="bar"){
    if(length(labels)>0){
      barplot(mat,names.arg=labels, cex.names = 0.9, xlab=xlab,ylab=ylab,main=main,ylim=range, col = colors, las=las)
    }else{
      barplot(mat,xlab=xlab,ylab=ylab,main=main, ylim=range, col=colors,cex.names = 0.9, las=las)
    }
  }else{
    stop(paste("Plot type",plotType,"is not supported."))
  }
  if(!is.null(colorMap) && (plotType!="pvalues") && (plotType!="violin")){
    legend("topright",legend=names(colorMap), cex=0.9, bg = "white", text.col=legend.colors)
  }
}


#' @title Assign colors such that they encode a metadata item
#' @description For a binary or categoric metadata item, as many colors as categories exist are assigned.
#' For a numeric metadata item, the item is binned according to the equalfreq or equalwidth
#' strategy in as many bins as requested. Alternatively, the user can also provide thresholds to define bins.
#' Each bin is assigned its color. If the groups vector is empty, a color vector is returned with as many entries as samples.
#' If groups are given, a group color vector is returned with as many entries as groups. The color map is returned along with
#' the color vector.
#' @details Notation: a round bracket means that the end point is excluded, a square bracket that the end point is included.
#' For groupColorStrategy maxfreq: if there are several equally frequent bins in the group, the first one is selected.
#' Missing values receive color "gray".
#' @param metadataItem a vector with binary, categoric or numeric values
#' @param groups an optional group vector with as many entries as samples in the target data set
#' @param binning binning strategy, only required for numeric metadata, either a string (equalfreq: the same number of samples in each group or equaldist: the same range of values) or a vector with thresholds
#' @param numBins the number of bins into which to bin a numeric metadata item, if zero: as many bins as samples (only required for numeric metadata and equalfreq/equaldist binning)
#' @param binLabels optional names for bins
#' @param groupColorStrategy maxfreq (assign the most frequent bin color to the group), binnum (color encodes number of bins per group)
#' @param returnBinned return the binned metadataItem
#' @return a list with color vector and a color map, if returnBinned is true, in addition bins
#' @export
makeColorsGivenMetadata<-function(metadataItem=c(), groups=c(), binning="equalfreq", numBins=0, binLabels=c(), groupColorStrategy="maxfreq",returnBinned=FALSE){

  bins=numBins
  if(bins==0){
    bins=length(metadataItem)
  }

  # vector with sample membership in binned metadata item
  binned=c()
  # vector with sample-wise or group-wise colors
  colors=c()
  # map to look up color assigned to bin
  colorMap=list()
  # result list
  res=list()
  # indices of missing values
  na.indices=which(is.na(metadataItem))

  # numeric metadata item and binning requested
  if(is.numeric(metadataItem) && bins<length(metadataItem)){
    print("Metadata item is numeric.")
    # binning strategy: thresholds were provided
    if(is.numeric(binning)){
      # check if minimum and maximum are provided as thresholds
      if(!(min(metadataItem) %in% binning)){
        binning=c(min(metadataItem),binning)
        print(paste("Adding minimum (",min(metadataItem),") as bin."))
      }
      if(!(max(metadataItem) %in% binning)){
        binning=c(binning,max(metadataItem))
        print(paste("Adding maximum (",max(metadataItem),") as bin."))
      }
      if(length(binLabels)>0){
        if(length(binLabels)==length(binning) || length(binLabels)>length(binning)){
          stop("There should be one label less than binning thresholds!")
        }
        binned=cut(metadataItem,breaks=binning, include.lowest=TRUE,labels=binLabels)
      }else{
        binned=cut(metadataItem,breaks=binning, include.lowest=TRUE)
      }
    # equaldist strategy
    }else if(binning == "equaldist"){
      # bin vector (equaldist)
      cut.out=cut(metadataItem,breaks=bins)
      print("Thresholds equaldist binning:")
      print(levels(cut.out))
      # rename bins (previously named after borders)
      if(length(binLabels)==bins){
        levels(cut.out)=binLabels
      }
      binned=cut.out
    # equalfreq strategy
    }else if(binning == "equalfreq"){
      # divide into requested number of intervals
      intervals=seq(0,1,(1/bins))
      thresholds=quantile(metadataItem,probs=intervals, na.rm=TRUE)
      print("Thresholds equalfreq binning:")
      print(thresholds[2:(length(thresholds)-1)])
      binned=metadataItem
      for(i in 2:length(thresholds)){
        set1=which(metadataItem>=thresholds[(i-1)])
        set2=which(metadataItem<thresholds[i])
        indices=intersect(set1,set2)
        if(length(binLabels)>0){
          binned[indices]=binLabels[(i-1)]
        }else{
          if(i==2){
            binned[indices]=paste("[",thresholds[(i-1)],",",thresholds[i],"]",sep="")
          }else if(i==length(thresholds)){
            binned[indices]=paste("[",thresholds[(i-1)],",",thresholds[i],"]",sep="")
          }else{
            binned[indices]=paste("(",thresholds[(i-1)],",",thresholds[i],")",sep="")
          }
        }
      }
      maxindex=which(metadataItem==max(metadataItem, na.rm=TRUE))
      if(length(binLabels)>0){
        binned[maxindex]=binLabels[length(binLabels)]
        levels(binned)=binLabels
      }else{
        binned[maxindex]=paste("[",thresholds[(length(thresholds)-1)],",",thresholds[length(thresholds)],"]",sep="")
      }
    }else{
      stop("Binning strategy ",binning," not supported!")
    }
  # do not bin numeric metadata item
  }else if(bins == length(metadataItem)){
    print("Each sample is assigned into its own bin.")
    binned = metadataItem
  }else{
    print("Metadata item is binary/categoric")
    bins=length(unique(metadataItem))
    print(paste("Number of categories:",bins))
    binned = metadataItem
  }

  # drop levels efficiently
  binned=as.character(binned)

  # treat missing values
  if(length(na.indices)>0){
    print(paste(length(na.indices),"missing values found."))
    binned[na.indices]="NA" # will be assigned a gray color
    colorMap[["NA"]]="gray"
  }

  #print(binned)

  # get a vector with as many colors as samples
  colors=assignColorsToGroups(groups=binned, refName="NA")

  # assemble color map
  for(index in 1:length(binned)){
    abin=binned[index]
    if(!(abin %in% names(colorMap))){
      colorMap[[as.character(abin)]]=colors[index]
    }
  }

  #print(names(colorMap))

  if(length(groups)>0){
    group.colors=c()
    unique.groups=unique(groups)
    hues=c()
    # count number of bins in each group to identify the range
    bins.per.group=c()
    group.counter=1
    if(groupColorStrategy=="binnum"){
      colorMap=list() # reset color map
      for(group in unique.groups){
        group.member.indices=which(groups==group)
        bins.per.group=c(bins.per.group,length(unique(binned[group.member.indices])))
      }
      col.vec = seq(0,1,1/max(bins.per.group))
      hues = hsv(col.vec)
    }
    # loop groups
    for(group in unique.groups){
      group.member.indices=which(groups==group)
      if(groupColorStrategy=="maxfreq"){
        # get the name of the most frequent bin in the group
        maxfreq.bin.label=names(sort(table(binned[group.member.indices]),decreasing = TRUE))[1]
        #print(paste("Most frequent bin in group",group,":",maxfreq.bin.label))
        # get an index of the bin
        maxfreq.bin.index=which(binned==maxfreq.bin.label)[1]
        # get the corresponding color
        maxfreq.bin.color=colors[maxfreq.bin.index]
        group.colors=c(group.colors,maxfreq.bin.color)
      }else if(groupColorStrategy=="binnum"){
        bin.number=bins.per.group[group.counter]
        group.colors=c(group.colors,hues[bin.number])
        if(!(bin.number %in% names(colorMap))){
          colorMap[[bin.number]]=hues[bin.number]
        }
      }else{
        stop(paste("Group color strategy",groupColorStrategy,"not supported."))
      }
      group.counter=group.counter+1
    } # end group loop
    colors=group.colors
  }

  if(returnBinned){
    res=list(colors, colorMap, binned)
    names(res)=c("colors","colormap","bins")
  }else{
    res=list(colors,colorMap)
    names(res)=c("colors","colormap")
  }
  return(res)
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
