#' @title Build a microbial network
#'
#' @description Wrapper for microbial network inference methods.
#' SPIEC-EASI and bnlearn are supported.
#' Abundances are assumed to have been already filtered and transformed.
#' The function returns the network as an igraph object,
#' with nodes labeled and colored at the selected taxonomic levels
#' and, for SPIEC-EASI, edges in red or green, depending on whether they
#' represent negative or positive associations. For SPIEC-EASI, edge width
#' represents the abolute of the regression coefficients, whereas for bnlearn, it represents
#' the probability of an arc to re-appear across bootstrap iterations.
#' If applicable, the color code for the node colors is plotted in a separate bar plot.
#' SPIEC-EASI returns undirected graphs, whereas bnlearn methods support the inference
#' of directed graphs. In general, SPIEC-EASI methods are faster than bnlearn methods.
#' The returned network can be visualized with igraph's interactive network visualisation
#' function tkplot or simply with igraph's plot function.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param lineages a matrix with lineages where row names match row names of abundances; first column gives kingdom, following columns give phylum, class, order, family, genus and species (obligatory for SPIEC-EASI)
#' @param method mb and glasso (SPIEC-EASI), hc, tabu, h2pc, mmhc, rsmax2, aracne, pc.stable, gs, fast.iamb, iamb.fdr, mmpc, hpc and si.hiton.pc (bnlearn, for aracne mi is set to mi-g)
#' @param repNum the number of bootstrap iterations in SPIEC-EASI and bnlearn; can be omitted (set to 0) for time-intensive bnlearn methods
#' @param minStrength bnlearn only: the minimum probability for an arc to appear across bootstraps (only applicable if repNum is set larger 1)
#' @param alpha p-value threshold for bnlearn (ignored in methods aracne, hc, tabu, rsmax2, mmhc and h2pc)
#' @param directed bnlearn only: infer a directed network (hc, tabu, rsmax2, mmhc and h2pc only infer directed networks, aracne only undirected networks)
#' @param nameLevel the taxonomic level used for node names (only if lineages provided)
#' @param colorLevel the taxonomic level used for node colors (only if lineages provided)
#' @param widthFactor edge width scaling factor
#' @return igraph object
#' @export

buildNetwork<-function(abundances,lineages=NULL,method="mb", repNum=20, minStrength=0.5, alpha=0.05, directed=FALSE, nameLevel="Genus", colorLevel="Class", widthFactor=10){

  supported.methods.bnlearn=c("hpc","si.hiton.pc","fast.iamb","gs","iamb.fdr","mmpc","pc.stable","hc","tabu","h2pc","mmhc","rsmax2","aracne")
  supported.methods.spiec=c("mb","glasso")
  result.graph=NULL

  otu.ids=rownames(abundances)

  # spiec-easi methods
  if(method %in% supported.methods.spiec){
    #strains=otu_table(abundances,taxa_are_rows = TRUE)
    #taxa=tax_table(lineages)
    #phyloseqobj=phyloseq(strains, taxa)
    phyloseqobj=toPhyloseq(abundances=abundances,lineages=lineages)
    spiec.out=spiec.easi(phyloseqobj, method=method,icov.select.params=list(rep.num=repNum))
    # get matrix of regression coefficients and symmetrize it
    betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))
    result.graph=graph_from_adjacency_matrix(betaMat,weighted=TRUE, mode="undirected")
    # does not allow building a weighted network, so is not used
    #result.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(phyloseqobj)))
    V(result.graph)$name=otu.ids
    edges=E(result.graph)
    edge.colors=c()
    edge.weights=c()
    print(paste("Number of edges found: ",length(edges)))
    if(length(edges)>0){
      # extract edge sign and assign edge color accordingly
      for(e.index in 1:length(edges)){
        adj.nodes=ends(result.graph,edges[e.index])
        xindex=which(otu.ids==adj.nodes[1])
        yindex=which(otu.ids==adj.nodes[2])
        beta=betaMat[xindex,yindex]
        edge.weights=c(edge.weights,abs(beta))
        if(beta>0){
          edge.colors=append(edge.colors,"forestgreen")
        }else if(beta<0){
          edge.colors=append(edge.colors,"red")
        }
      }
      E(result.graph)$color=edge.colors
      E(result.graph)$weight=edge.weights
      E(result.graph)$width=edge.weights*widthFactor
      result.graph=colorGraph(result.graph,otu.ids = otu.ids,lineages = lineages,colorLevel = colorLevel,nameLevel = nameLevel)
      # remove orphan nodes (not necessary for bnlearn)
      result.graph=delete.vertices(result.graph,degree(result.graph)==0)
    }
  }else{
    bnlearn.data=as.data.frame(t(abundances))
    # argument B not required in any of the methods listed here
    if(method=="hpc"){
      bnlearn.out=hpc(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="gs"){
      bnlearn.out=gs(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="si.hiton.pc"){
      bnlearn.out=si.hiton.pc(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="iamb.fdr"){
      bnlearn.out=iamb.fdr(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="pc.stable"){
      bnlearn.out=pc.stable(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="fast.iamb"){
      bnlearn.out=fast.iamb(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="mmpc"){
      bnlearn.out=mmpc(bnlearn.data, alpha=alpha, undirected = (!directed))
    }else if(method=="aracne"){
      directed=FALSE
      bnlearn.out=aracne(bnlearn.data,mi="mi-g")
    }else if(method=="hc"){
      directed=TRUE
      bnlearn.out=hc(bnlearn.data)
    }else if(method=="tabu"){
      directed=TRUE
      bnlearn.out=tabu(bnlearn.data)
    }else if(method=="h2pc"){
      directed=TRUE
      bnlearn.out=h2pc(bnlearn.data)
    }else if(method=="rsmax2"){
      bnlearn.out=rsmax2(bnlearn.data)
    }else if(method=="mmhc"){
      directed=TRUE
      bnlearn.out=mmhc(bnlearn.data)
    }else{
      stop(paste("Method not supported. Please choose a bnlearn method (",paste(supported.methods.bnlearn,collapse=", "),") or a SPIEC-EASI method (",paste(supported.methods.spiec,collapse=", "),").",sep=""))
    }
    print(bnlearn.out)

    if(repNum>1){
      # arc.strength does not require bootstrapping, but it only works for entirely directed methods
      print(paste("Carrying out ",repNum,"bootstraps to assess edge strengths..."))
      arcs=boot.strength(bnlearn.data,R=repNum,algorithm=method)
      edge.weights=as.numeric(arcs[,3]) # strength
      indices.selected.edges=which(edge.weights>minStrength)
      #print(indices.selected.edges)
      if(length(indices.selected.edges)==0){
        stop("All edges have a strength below the selected minimum strength.")
      }
      result.graph=graph_from_edgelist(as.matrix(arcs[indices.selected.edges,1:2]), directed=directed)
      E(result.graph)$weight=edge.weights[indices.selected.edges]
      E(result.graph)$width=edge.weights[indices.selected.edges]*widthFactor
    }else{
      result.graph=graph_from_edgelist(bnlearn.out$arcs, directed=directed)
      warning("Please select bootstrap iterations to set edge strengths with bnlearn.")
    }

    # if graph is entirely directed, report overall goodness of fit and extract edge strengths
    if(directed(bnlearn.out)){
      fitted=bn.fit(bnlearn.out,bnlearn.data)
      print(paste("logLik:",logLik(fitted,bnlearn.data)))
      print(paste("AIC:",AIC(fitted,bnlearn.data)))
      print(paste("BIC:",BIC(fitted,bnlearn.data)))
    }

    result.graph=colorGraph(result.graph,otu.ids = otu.ids,lineages = lineages,colorLevel = colorLevel,nameLevel = nameLevel)
  }

  return(result.graph)
}

# label and color nodes in igraph object according to selected lineage level and set global attributes
colorGraph<-function(result.graph, otu.ids=c(),lineages=NULL,nameLevel="",colorLevel=""){
  # lineages provided
  if(!is.null(lineages)){
    nodenames=getTaxonomy(otu.ids, lineages, useRownames=TRUE, level=nameLevel)
    # nodes belonging to the same taxonomic level receive the same node color
    colorgroups=getTaxonomy(otu.ids, lineages, useRownames=TRUE, level=colorLevel)
    treated.colorgroups=c()
    # treat missing taxonomic assignments
    for(color.member in colorgroups){
      if(is.na(color.member)){
        color.member="NA"
      }
      treated.colorgroups=c(treated.colorgroups,color.member)
    }
    colorgroups=treated.colorgroups
    assignedColors=assignColorsToGroups(colorgroups, returnMap=TRUE)
    nodecolors=assignedColors$colors
    #print(assignedColors$colorMap)
    V(result.graph)$name=nodenames
    V(result.graph)$color=nodecolors
    unique.colors=c()
    unique.names=c()
    for(name in names(assignedColors$colorMap)){
      unique.colors=c(unique.colors,assignedColors$colorMap[[name]])
      unique.names=c(unique.names,name)
    }
    values=rep(1,length(unique.names))
    pre.par=par() # save par values
    par(mar=c(2,14,2,2)) # adjust margins for barplot
    # plot node color legend
    barplot(values, horiz=TRUE, names.arg = unique.names,las=2,xlim=c(0,10), xaxt="n", col=unique.colors, main=paste(colorLevel,"legend"))
    par=pre.par # restore default par
  }else{
    V(result.graph)$color="white"
  }
  E(result.graph)$arrow.size=1 # previously 5, but too big now
  V(result.graph)$vertex.size=5
  V(result.graph)$frame.color="black"
  V(result.graph)$label.color="black"
  return(result.graph)
}



#' @title Get the taxonomy given OTU names and lineage information
#'
#' @description The lineage information is provided in form of a matrix,
#' which contains for each OTU identifier the taxonomic levels from column
#' 2 to 8 (domain, phylum, class, order, family, genus, species) and the entire
#' lineage in column 9. The lineage in column 9 is optional.
#' Alternatively, OTU identifiers can also be stored as column names.
#' In this case, taxonomic levels are assumed to range from column 1 to 7.
#' To differentiate between these two formats, set useRownames to false or true.
#'
#' @param selected a character vector of selected OTU identifiers
#' @param lineages a lineage table
#' @param level the taxonomic level (domain, phylum, class, order, family, genus or species)
#' @param useRownames match OTU identifiers to row names instead of the first column (in this case, taxonomic levels are assumed to range from column 1 to 7)
#' @return the taxonomy of the OTUs
getTaxonomy<-function(selected=c(),lineages,level="class", useRownames=FALSE){
  # exact OTU identifier match
  taxa=c()
  domainIndex=2
  if(useRownames==TRUE){
    domainIndex=1
    indices=match(selected,rownames(lineages))
  }else{
    indices=match(selected,lineages[,1])
  }
  if(tolower(level)=="domain"){
    taxa=lineages[indices,domainIndex]
  }else if(tolower(level)=="phylum"){
    taxa=lineages[indices,(domainIndex+1)]
  }else if(tolower(level)=="class"){
    taxa=lineages[indices,(domainIndex+2)]
  }
  else if(tolower(level)=="order"){
    taxa=lineages[indices,(domainIndex+3)]
  }
  else if(tolower(level)=="family"){
    taxa=lineages[indices,(domainIndex+4)]
  }
  else if(tolower(level)=="genus"){
    taxa=lineages[indices,(domainIndex+5)]
  }
  else if(tolower(level)=="species"){
    taxa=lineages[indices,(domainIndex+6)]
  }
  return(taxa)
}
