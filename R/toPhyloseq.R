#' @title Convert abundances to phyloseq object
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param metadata an optional data frame with metadata items as columns, where samples are in the same order as in abundances
#' @param lineages an optional lineage matrix where row names match row names of taxa; taxa must be provided in the same order as in x; first column gives kingdom, following columns give phylum, class, order, family, genus and species
#' @param assignSampleNames assign random sample names
#' @param dummyMetadata create a dummy metadata object with a single numeric item from a normal distribution; only if no metadata are provided
#' @param dummyLineages assign dummy lineages; only if no lineages are provided and a lineageSource is given
#' @param lineageSource a phyloseq object with known lineages, e.g. GlobalPatterns (phyloseq preloaded data that can be loaded with data(GlobalPatterns)), to be used for random lineage assignment
#' @return phyloseq object
#' @examples
#' data(ibd_taxa)
#' data(ibd_metadata)
#' data(ibd_lineages)
#' phyloseq.obj=toPhyloseq(ibd_taxa,metadata=ibd_metadata,lineages=ibd_lineages)
#' phyloseq.obj
#' @export
toPhyloseq<-function(abundances, metadata=NULL, lineages=NULL, assignSampleNames=FALSE, dummyMetadata=FALSE, dummyLineages=FALSE, lineageSource=NULL){

  if(!(is.null(metadata))){
    if(ncol(abundances)!=nrow(metadata)){
     stop("There should be as many metadata rows as there are columns in abundances!")
    }
    if(is.matrix(metadata)){
      metadata=as.data.frame(metadata)
    }
  }

  if(!is.null(lineages)){
    lineages=as.matrix(lineages)
  }

  # assign taxon names if necessary
  if(is.null(rownames(abundances))){
    otus=c()
    for(i in 1:nrow(abundances)){
      otus=c(otus,paste("OTU",i,sep=""))
    }
    rownames(abundances)=otus
  }

  # assign random sample names if necessary
  if(is.null(colnames(abundances)) || assignSampleNames){
    colnames.rand=c()
    for(i in 1:ncol(abundances)){
      colnames.rand=c(colnames.rand,makeRandomString())
    }
    colnames(abundances)=colnames.rand
  }

  # instantiate phyloseq objects
  phylo_taxa=otu_table(abundances,taxa_are_rows = TRUE)
  phylo_sample_data=NULL
  phylo_lineages=NULL

  # create random sample data
  if(is.null(metadata) && dummyMetadata){
    random.metadata=rnorm(ncol(abundances))
    sam.dummy=as.matrix(random.metadata)
    rownames(sam.dummy)=colnames(abundances)
    colnames(sam.dummy)=c("dummy")
    phylo_sample_data=sample_data(as.data.frame(sam.dummy))
  }
  # convert provided metadata
  if(!is.null(metadata)){
    phylo_sample_data=sample_data(metadata)
  }

  # create random lineages
  if(is.null(lineages) && dummyLineages){
    if(is.null(lineageSource)){
      stop("Please provide a phyloseq object with lineages, e.g. data(GlobalPatterns).")
    }else{
      lineage.mat=matrix(nrow=nrow(abundances),ncol=7)
      temp.lineages=tax_table(lineageSource)
      for(i in 1:nrow(abundances)){
          lineage.mat[i,]=getRandomLineage(temp.lineages)
      }
      rownames(lineage.mat)=rownames(abundances)
      phylo_lineages=tax_table(lineage.mat)
    }
  }
  # convert provided lineages
  if(!is.null(lineages)){
    if(ncol(lineages)!=7){
      stop("Lineage matrix needs to have 7 columns for kingdom, phylium, class, order, family, genus and species.")
    }
    if(nrow(lineages)!=nrow(abundances)){
      stop("Row number of lineage matrix differ from row number in abundance matrix!")
    }
    if(is.null(rownames(lineages))){
      warning("Lineage matrix does not provide row names. Row names from abundance matrix are assigned.")
      rownames(lineages)=rownames(abundances)
    }
    phylo_lineages=tax_table(lineages)
  }

  # return phyloseq object
  if(is.null(phylo_sample_data) && is.null(phylo_lineages)){
    return(phyloseq(phylo_taxa))
  }else if(is.null(phylo_sample_data)){
    return(phyloseq(phylo_taxa,phylo_lineages))
  }else if(is.null(phylo_lineages)){
    return(phyloseq(phylo_taxa,phylo_sample_data))
  }else{
    return(phyloseq(phylo_taxa,phylo_sample_data, phylo_lineages))
  }
}

# get a random entry from the tempLineage tax_table
getRandomLineage<-function(tempLineages){
  index=sample(c(1:nrow(tempLineages)))[1]
  return(tempLineages[index,])
}

# Function to generate random strings adapted from:
# https://ryouready.wordpress.com/2008/12/18/generate-random-string-name/
makeRandomString <- function(length=10){
    return(paste(sample(c(0:9, letters, LETTERS),length, replace=TRUE),collapse=""))
}
