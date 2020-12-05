#' @title Build a microbial network
#'
#' @description Wrapper for microbial network inference methods.
#' SPIEC-EASI, CoNet (see function \code{\link{barebonesCoNet}}) and bnlearn are supported.
#' Abundances are assumed to have been already filtered and transformed.
#' The function returns the network as an igraph object,
#' with nodes labeled and colored at the selected taxonomic levels
#' and, for SPIEC-EASI and CoNet, edges in red or green, depending on whether they
#' represent negative or positive associations.
#'
#' @details For SPIEC-EASI, edge width represents the abolute of the regression coefficients, whereas for bnlearn,
#' it represents the probability of an arc to re-appear across bootstrap iterations. For CoNet, it either
#' represents absolute association strength, method number supporting an edge or significance,
#' depending on whether more than one method was used to compute the network and whether or not
#' p-values were computed (significance is defined as -log10(pval)). For single dissimilarities in CoNet, edge weight
#' is scaled (see \code{\link{barebonesCoNet}} for details) such that for all network construction methods,
#' a thicker edge means a stronger/more stable/more significant association.
#' For CoNet, p-values are computed by default with cor.test for correlations and permutations and bootstraps
#' for dissimilarities and are multiple-testing corrected with Benjamini-Hochberg.
#' P-value computation in CoNet can be disabled by setting repNum to 0.
#' If lineages are provided and legend is true, the color code for the node colors is plotted in a separate bar plot.
#' SPIEC-EASI and CoNet return undirected graphs, whereas bnlearn methods support the inference
#' of directed graphs. SPIEC-EASI methods tend to be faster than bnlearn methods.
#' The returned network can be visualized with igraph's interactive network visualisation
#' function tkplot or simply with igraph's plot function. When RCy3 is installed and Cytoscape is open,
#' the graph can also be sent to Cytoscape directly with RCy3's function createNetworkFromIgraph().
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param lineages a matrix with lineages where row names match row names of abundances; first column gives kingdom, following columns give phylum, class, order, family, genus and species (obligatory for SPIEC-EASI)
#' @param min.occ only keep rows with at least the given number of non-zero values, keep the sum of the other rows (carried out before network construction)
#' @param norm normalize matrix (carrried out after filtering); not carried out for SPIEC-EASI (where it is already applied)
#' @param clr apply CLR transform (after filtering and normalization; \code{\link{clr}} with omit.zeros true); not carried out for SPIEC-EASI (already applied)
#' @param method mb or glasso (SPIEC-EASI), any combination of pearson, spearman, bray and/or kld (CoNet), hc, tabu, h2pc, mmhc, rsmax2, aracne, pc.stable, gs, fast.iamb, iamb.fdr, mmpc, hpc or si.hiton.pc (bnlearn, for aracne mi is set to mi-g)
#' @param repNum the number of permutation/bootstrap iterations for CoNet with dissimilarities and of bootstrap iterations in SPIEC-EASI and bnlearn; can be omitted (set to 0) for time-intensive bnlearn methods or for CoNet without p-values
#' @param minStrength bnlearn only: the minimum probability for an arc to appear across bootstraps (only applicable if repNum is set larger 1)
#' @param alpha p-value threshold for bnlearn and CoNet (ignored in methods aracne, hc, tabu, rsmax2, mmhc and h2pc)
#' @param renorm compute correlation p-values in CoNet with renorm enabled (slow)
#' @param initEdgeNum CoNet parameter for the number of initial top and bottom edges
#' @param methodNum CoNet threshold on the minimum method number needed to keep an edge (only used when more than one method is provided)
#' @param directed bnlearn only: infer a directed network (hc, tabu, rsmax2, mmhc and h2pc only infer directed networks, aracne only undirected networks)
#' @param nameLevel the taxonomic level used for node names (only if lineages provided)
#' @param colorLevel the taxonomic level used for node colors (only if lineages provided)
#' @param widthFactor edge width scaling factor
#' @param legend if lineages are provided: show node color code
#' @param left.margin if lineages are provided and legend is true: set the left margin (the margin is restored to its previous value afterwards)
#' @return igraph object
#' @examples
#' data(ibd_taxa)
#' data(ibd_lineages)
#' top.indices=sort(rowSums(ibd_taxa),decreasing = TRUE,index.return=TRUE)$ix[1:30]
#' taxa=ibd_taxa[top.indices,]
#' lineages=ibd_lineages[top.indices,]
#' hc.out=buildNetwork(taxa,lineages,method="hc",nameLevel="species")
#' plot(hc.out)
#' @export

buildNetwork<-function(abundances, lineages=NULL, min.occ=0, norm=FALSE, clr=FALSE, method="mb", repNum=20, minStrength=0.5, alpha=0.05, renorm=FALSE, initEdgeNum=100, methodNum=2, directed=FALSE, nameLevel="Genus", colorLevel="Class", widthFactor=10, legend=TRUE, left.margin=10){

  supported.methods.bnlearn=c("hpc","si.hiton.pc","fast.iamb","gs","iamb.fdr","mmpc","pc.stable","hc","tabu","h2pc","mmhc","rsmax2","aracne")
  supported.methods.spiec=c("mb","glasso")
  supported.methods.conet=c("pearson","spearman","kld","bray")
  result.graph=NULL

  if(length(method)==1 && method %in% supported.methods.spiec && is.null(lineages)){
    stop("For SPIEC-EASI, lineages are obligatory!")
  }

  # CoNet methods can be longer 1
  if(length(method)==1){
    if(method %in% supported.methods.spiec && (norm==TRUE || clr==TRUE)){
      norm=FALSE
      clr=FALSE
      warning("SPIEC-EASI normalizes and CLR-transforms data, so norm and clr are ignored.")
    }
  }

  if(min.occ>0){
    abundances = filterTaxonMatrix(abundances,keepSum = TRUE, minocc=min.occ)
    # filter lineages if given
    if(!is.null(lineages)){
      indices.kept=match(rownames(abundances),rownames(lineages))
      lineages=lineages[indices.kept,]
    }
  }

  # normalize matrix
  if(norm == TRUE){
    abundances=normalize(abundances)
  }

  if(clr){
    # filtering has been carried out before
    abundances=clr(abundances = abundances,minocc = 0, omit.zeros = TRUE)
  }

  otu.ids=rownames(abundances)

  # spiec-easi methods
  if(length(method)==1 && method %in% supported.methods.spiec){
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
      # remove orphan nodes (not necessary for bnlearn)
      result.graph=delete.vertices(result.graph,degree(result.graph)==0)
    }
  }else if(length(method)>1 || method %in% supported.methods.conet){
      pval.cor=FALSE
      permut=TRUE
      permutandboot=TRUE
      if(repNum==0 || is.na(repNum)){
        permut=FALSE
        pval.cor=FALSE
        permutandboot=FALSE
      }
      if("spearman" %in% method || "pearson" %in% method){
        if(!renorm && repNum>0){
          pval.cor=TRUE
        }
      }
      result.graph=barebonesCoNet(abundances = abundances, methods = method, pval.T = alpha, bh=TRUE, min.occ = min.occ, init.edge.num = initEdgeNum, method.num.T = methodNum,iters = repNum, renorm=renorm, permut=permut,permutandboot = permutandboot, pval.cor = pval.cor)
      #print("graph built")
      E(result.graph)$width=E(result.graph)$weight*widthFactor
      # removal of orphan nodes not necessary for CoNet
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
      stop(paste("Method not supported. Please choose a bnlearn method (",paste(supported.methods.bnlearn,collapse=", "),"), a CoNet method (",paste(supported.methods.conet,collapse=", "),") or a SPIEC-EASI method (",paste(supported.methods.spiec,collapse=", "),").",sep=""))
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
    if(bnlearn::directed(bnlearn.out)){
      fitted=bn.fit(bnlearn.out,bnlearn.data)
      print(paste("logLik:",logLik(fitted,bnlearn.data)))
      print(paste("AIC:",AIC(fitted,bnlearn.data)))
      print(paste("BIC:",BIC(fitted,bnlearn.data)))
    }
  }

  result.graph=colorGraph(result.graph,lineages = lineages,colorLevel = colorLevel,nameLevel = nameLevel, legend=legend, left.margin=left.margin)
  return(result.graph)
}

# label and color nodes in igraph object according to selected lineage level and set global attributes
colorGraph<-function(result.graph, lineages=NULL,nameLevel="",colorLevel="", legend=FALSE, left.margin=10){
  # if graph is non-empty
  if(length(E(result.graph))>0){
    # lineages provided
    if(!is.null(lineages)){
      otu.ids=V(result.graph)$name
      #print(otu.ids[1])
      #print(otu.ids[1:10])
      #print(nameLevel)
      #print(colorLevel)
      #print(lineages[1,])
      nodenames=getTaxonomy(otu.ids, lineages, useRownames=TRUE, level=nameLevel)
      # nodes belonging to the same taxonomic level receive the same node color
      colorgroups=getTaxonomy(otu.ids, lineages, useRownames=TRUE, level=colorLevel)
      #print(colorgroups)
      treated.colorgroups=c()
      # treat missing taxonomic assignments
      for(color.member in colorgroups){
        if(is.na(color.member)){
          color.member="NA"
        }
        treated.colorgroups=c(treated.colorgroups,color.member)
      }
      colorgroups=treated.colorgroups
      #print(length(colorgroups))
      assignedColors=assignColorsToGroups(colorgroups, returnMap=TRUE)
      nodecolors=assignedColors$colors
      #print(assignedColors$colorMap)
      #print(length(V(result.graph)))
      #print(length(nodenames))
      V(result.graph)$name=nodenames
      #print(length(nodecolors))
      #print(length(V(result.graph)))
      V(result.graph)$color=nodecolors
      unique.colors=c()
      unique.names=c()
      for(name in names(assignedColors$colorMap)){
        unique.colors=c(unique.colors,assignedColors$colorMap[[name]])
        unique.names=c(unique.names,name)
      }
      values=rep(1,length(unique.names))
      if(legend){
        pre.mar=par()$mar
        # bottom, left, top, and right
        par(mar=c(2,left.margin,2,2)) # adjust margins for barplot
        # plot node color legend
        barplot(values, horiz=TRUE, names.arg = unique.names,las=2,xlim=c(0,10), xaxt="n", col=unique.colors, main=paste(colorLevel,"legend"))
        par(mar=pre.mar)
      }
    }else{
      V(result.graph)$color="white"
    }
    E(result.graph)$arrow.size=1 # previously 5, but this is too big
    V(result.graph)$vertex.size=5
    V(result.graph)$frame.color="black"
    V(result.graph)$label.color="black"
  }
  return(result.graph)
}

# Get the taxonomy given OTU names and lineage information
#
# The lineage information is provided in form of a matrix,
# which contains for each OTU identifier the taxonomic levels from column
# 2 to 8 (domain, phylum, class, order, family, genus, species) and the entire
# lineage in column 9. The lineage in column 9 is optional.
# Alternatively, OTU identifiers can also be stored as column names.
# In this case, taxonomic levels are assumed to range from column 1 to 7.
# To differentiate between these two formats, set useRownames to false or true.
#
# selected a character vector of selected OTU identifiers
# lineages a lineage table
# level the taxonomic level (domain, phylum, class, order, family, genus or species)
# useRownames match OTU identifiers to row names instead of the first column (in this case, taxonomic levels are assumed to range from column 1 to 7)
# returns the taxonomy of the OTUs
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
