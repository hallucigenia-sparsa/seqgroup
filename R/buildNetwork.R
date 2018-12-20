#' @title Build a microbial network
#'
#' @description Wrapper for microbial network inference methods.
#' For the moment, only SPIEC-EASI is supported.
#' Abundances are assumed to have been already filtered.
#' The function returns the network as an igraph object,
#' with nodes labeled and colored at the selected taxonomic levels
#' and edges in red or green, depending on whether they
#' represent negative or positive associations. Edge thickness
#' represents association strength. The color code for the node colors is
#' plotted in a separate bar plot.
#' The returned network can be visualized with igraph's
#' interactive network visualisation function tkplot or
#' simply with igraph's plot function.
#' TODO: support bnlearn
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param lineages a matrix with lineages in the phyloseq format
#' @param method mb and glasso, supported by SPIEC-EASI
#' @param repNum the number of bootstrap iterations in SPIEC-EASI
#' @param nameLevel the taxonomic level used for node names
#' @param colorLevel the taxonomic level used for node colors
#' @param widthFactor edge width scaling factor (for SPIEC-EASI, the width is the abolute of the regression coefficients)
#' @return the SPIEC-EASI igraph object
#' @export

buildNetwork<-function(abundances,lineages,method="mb", repNum=20, nameLevel="Genus", colorLevel="Class", widthFactor=50){
  strains=otu_table(abundances,taxa_are_rows = TRUE)
  taxa=tax_table(lineages)
  phyloseqobj=phyloseq(strains, taxa)
  spiec.out=spiec.easi(phyloseqobj, method=method,icov.select.params=list(rep.num=repNum))
  # get matrix of regression coefficients and symmetrize it
  betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))
  spiec.graph=graph_from_adjacency_matrix(betaMat,weighted=TRUE, mode="undirected")
  # does not allow building a weighted network, so is not used
  #spiec.graph=adj2igraph(getRefit(spiec.out), vertex.attr=list(name=taxa_names(phyloseqobj)))
  otu.ids=rownames(strains)
  V(spiec.graph)$name=otu.ids
  edges=E(spiec.graph)
  edge.colors=c()
  edge.weights=c()
  print(paste("Number of edges found: ",length(edges)))
  if(length(edges)>0){
    # extract edge sign and assign edge color accordingly
    for(e.index in 1:length(edges)){
      adj.nodes=ends(spiec.graph,edges[e.index])
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
    E(spiec.graph)$color=edge.colors
    E(spiec.graph)$weight=edge.weights
    E(spiec.graph)$width=edge.weights*widthFactor
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
    V(spiec.graph)$name=nodenames
    E(spiec.graph)$arrow.size=5
    V(spiec.graph)$color=nodecolors
    V(spiec.graph)$frame.color="black"
    V(spiec.graph)$label.color="black"
    unique.colors=c()
    unique.names=c()
    for(name in names(assignedColors$colorMap)){
      unique.colors=c(unique.colors,assignedColors$colorMap[[name]])
      unique.names=c(unique.names,name)
    }
    values=rep(1,length(unique.names))
    pre.par=par()
    par(mar=c(2,14,2,2))
    # plot node color legend
    barplot(values, horiz=TRUE, names.arg = unique.names,las=2,xlim=c(0,10), xaxt="n", col=unique.colors, main=paste(colorLevel,"legend"))
    par=pre.par
    # remove orphan nodes
    spiec.graph=delete.vertices(spiec.graph,degree(spiec.graph)==0)
  }
  return(spiec.graph)
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
