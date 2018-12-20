#' @title Export taxon abundances in a format suitable for CoNet
#'
#' @description If provided, metadata are exported into a CoNet feature file. Prior to the export, categoric
#' metadata are binarised, such that the resulting feature file only contains binary (as 0/1) and
#' numeric features, since CoNet cannot handle categoric features. Constant metadata are removed.
#' If binary metadata are provided as strings, the strings encoding 0 (no) and 1 (yes) have to be indicated.
#' If lineages are provided, they are exported into a CoNet metadata file, which can be used for all group-specific
#' data sets, in case groups are provided. Either the last column or the row names of the lineage matrix have to match the
#' row names of the abundance matrix.
#'
#' TODO: treat special characters
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param metadata an optional data frame with metadata items where sample names match sample names in abundances; data types (factor vs numeric) are supposed to be correct
#' @param lineages an optional lineage matrix where row names match row names of taxa; taxa must be provided in the same order as in abundances; first column gives kingdom, following columns give phylum, class, order, family, genus and species (not all are required)
#' @param groups an optional group membership vector with as many entries as abundances has samples; if provided, data are exported separately for each group
#' @param export.folder the export folder into which files should be saved
#' @param root.name a name included in all exported files
#' @param metadata.to.skip names of metadata items that should not be exported
#' @param omitNeg if true, metadata with negative values are removed, since CoNet run with Bray Curtis or Kullback Leibler dissimilarity cannot handle them
#' @param yes the string encoding 1 in a binary metadata item
#' @param no the string encoding 0 in a binary metadata item
#' @param date.items names of metadata items that are dates (they will be converted into days since 1/1/1900)
#' @param date.format the date format (an example date fitting the default format is 26/1/80)
#' @param taxa.are.rownames if true, row names instead of the last column of the lineage matrix provide the lowest-level taxa (e.g. for OTUs or sequencing variants)
#' @export

exportToCoNet<-function(abundances, metadata=NULL, lineages=NULL, groups=c(), export.folder="", root.name="conet", metadata.to.skip=c(), omitNeg=TRUE, yes="Y", no="N", date.items=c(), date.format="%d/%m/%y", taxa.are.rownames=FALSE){
  if(is.data.frame(metadata)==FALSE){
    metadata=as.data.frame(metadata)
  }
  for(skip in metadata.to.skip){
    metadata[[skip]]=NULL
  }
  if(length(groups)>0){
    unique.groups=unique(groups)
    for(group in unique.groups){
    export.abundance.file=file.path(export.folder,paste(root.name,"_",group,"_abundances.txt",sep=""))
    export.metadata.file=file.path(export.folder,paste(root.name,"_",group,"_features.txt",sep=""))
    group.indices=which(groups==group)
    group.abundances=abundances[,group.indices]
    group.metadata=metadata
    if(!is.null(metadata)){
      group.metadata=metadata[group.indices,] # in the data frame, metadata are columns and samples rows
    }
    exportSubset(abundances=group.abundances,metadata=group.metadata,export.abundance.file = export.abundance.file,export.metadata.file = export.metadata.file, omitNeg=omitNeg, yes=yes, no=no, date.items=date.items, date.format=date.format)
    }
  }else{
    export.abundance.file=file.path(export.folder,paste(root.name,"_abundances.txt",sep=""))
    export.metadata.file=file.path(export.folder,paste(root.name,"_features.txt",sep=""))
    exportSubset(abundances=abundances,metadata=metadata,export.abundance.file = export.abundance.file,export.metadata.file = export.metadata.file, omitNeg=omitNeg, yes=yes, no=no, date.items=date.items, date.format=date.format)
  }
  # export lineages as metadata file (it is not group-specific)
  if(!is.null(lineages)){
    lineages.export.table=formatLineages(lineages, taxa.in.rownames = taxa.are.rownames)
    export.lineage.file=file.path(export.folder,paste(root.name,"_metadata.txt",sep=""))
    write.table(lineages.export.table, file=export.lineage.file, row.names=FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
}

# export abundances and, if provided, metadata
exportSubset<-function(abundances, metadata=NULL, export.abundance.file="", export.metadata.file="", omitNeg=TRUE, yes="Y", no="N", date.items=c(), date.format="%d/%m/%y"){
  # CoNet requires a name for the row name column
  row.names=rownames(abundances) # taxon names
  col.names=colnames(abundances) # sample names
  abundances=cbind(row.names,abundances)
  colnames(abundances)=c("Taxa", col.names)
  write.table(abundances, file=export.abundance.file, row.names=FALSE, quote = FALSE, sep = "\t", na = "NaN")
  # export metadata
  if(!is.null(metadata)){
    # make metadata numeric (also removes constant and negative metadata)
    # note that some binary data may have white spaces, which will be removed by CoNet
    num.metadata=metadataToNumeric(metadata,yes=yes,no=no, date.items = date.items, format=date.format, remove.neg=omitNeg)
    # introduce a metadata item with metadata names
    metadata.names=names(num.metadata) # metadata names
    num.metadata=rbind(NAME=metadata.names,num.metadata)
    write.table(t(num.metadata), file=export.metadata.file, row.names=FALSE, quote = FALSE, sep = "\t", na = "NaN")
  }
}

# reformat lineages as a CoNet metadata table
# If taxa.in.rownames is true, the last column in the lineage file does not
# correspond to the lowest taxonomic level, which is instead encoded in the row names
# this is the case for instance for OTUs and sequencing variants
# missing lineage information may be encoded as NA
formatLineages<-function(lineages, taxa.in.rownames=FALSE){
  lineages.export.table=matrix(nrow=nrow(lineages),ncol=3)
  levelNum=ncol(lineages)
  for(row.index in 1:nrow(lineages)){
    if(taxa.in.rownames){
      taxon=rownames(lineages)[row.index]
    }else{
      taxon=lineages[row.index,levelNum]
    }
    # start with kingdom
    lineage=lineages[row.index,1]
    levelCounter=2
    # append all other levels
    while(!is.na(lineages[row.index,levelCounter]) && levelCounter<levelNum){
      lineage=paste(lineage,lineages[row.index,levelCounter],sep="--")
      levelCounter=levelCounter+1
    }
    lineage=paste(lineage,taxon,sep="--")
    #print(lineage)
    lineages.export.table[row.index,1]=taxon
    lineages.export.table[row.index,2]=lineage
    lineages.export.table[row.index,3]=taxon
  }
  return(lineages.export.table)
}
