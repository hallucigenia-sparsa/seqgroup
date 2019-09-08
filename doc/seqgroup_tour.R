## ----setup, include=FALSE------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_seqgroup_tour/")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(seqgroup)

## ---- message=FALSE, warning=FALSE---------------------------------------
data("ibd_taxa")
data("ibd_lineages")
data("ibd_metadata")
summary(ibd_metadata)

## ---- message=FALSE, warning=FALSE---------------------------------------
groups=as.vector(ibd_metadata$Diagnosis)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
compareGroups(ibd_taxa,groups=groups,property="alpha",pvalViz = TRUE)
compareGroups(ibd_taxa,groups=groups,property="beta",pvalViz = TRUE)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 3, fig.width = 3----
compareGroups(ibd_taxa,groups=groups,property="beta",method="DM")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
control.indices=which(groups=="Control")
groupBarplot(ibd_taxa[,control.indices],topTaxa = 15,legend.shift = 1.5, main="Control")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
groupBarplot(ibd_taxa,topTaxa = 15,groups=groups, randSampleNum = 10, legend.shift = 1.5, main="Taxon composition")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
ibd.indices=c(which(groups=="UC"),which(groups=="CD"))
par(mfrow=c(2,1))
taxon.color.map=groupBarplot(ibd_taxa[,control.indices],topTaxa = 7, randSampleNum = 20, extendTaxonColorMap = TRUE, main="Control")
taxon.color.map=groupBarplot(ibd_taxa[,ibd.indices],topTaxa = 7, randSampleNum = 20, taxon.color.map = taxon.color.map, extendTaxonColorMap = TRUE, main="IBD")
par(mfrow=c(1,1))

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
library(phyloseq)
prism.obj=toPhyloseq(ibd_taxa,metadata=ibd_metadata,lineages=ibd_lineages)
p = plot_bar(prism.obj, "Order", fill="Order", facet_grid=~Diagnosis) + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")
p

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
seqPCoA(ibd_taxa,groups=groups, topTaxa=30)
# prevalence filter
prev=50 # minimum occurrence in 50 samples
ibd_taxa_incidence=ibd_taxa
ibd_taxa_incidence[ibd_taxa_incidence>0]=1
rowSums=rowSums(ibd_taxa_incidence)
indices.prev=which(rowSums>=prev)
seqPCoA(ibd_taxa[indices.prev,], groups=groups, topTaxa=30)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
# assign numeric metadata types, metadata with only 2 values are treated automatically as categoric (binary).
ibd_metadata=assignMetadataTypes(ibd_metadata,categoric=c("SRA_metagenome_name","Diagnosis"))
# remove 65 samples with missing calprotectin measurements or other missing values in the metadata
na.indices=unique(which(is.na(ibd_metadata),arr.ind=TRUE)[,1])
indices.to.keep=setdiff(1:nrow(ibd_metadata),na.indices)
ibd_metadata.nafree=ibd_metadata[indices.to.keep,]
ibd_taxa.nafree=ibd_taxa[,indices.to.keep]
seqPCoA(ibd_taxa.nafree,metadata=ibd_metadata.nafree,groups=groups[indices.to.keep],topTaxa=60)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
seqPCoA(ibd_taxa.nafree,metadata=ibd_metadata.nafree,groups=groups[indices.to.keep],topTaxa=60, dis="gower",xlim=c(-0.3,0.3),main="Gower")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
ibd.groups=groups[indices.to.keep]
ibd.groups[ibd.groups=="UC"]="IBD"
ibd.groups[ibd.groups=="CD"]="IBD"
taxon="Faecalibacterium_prausnitzii"
compareDistribs(ibd_taxa.nafree,taxon=taxon,groups=ibd.groups,group1="IBD",group2="Control")

## ---- message=FALSE, results='hide', warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
clusters=findClusters(ibd_taxa,k=NA,maxK=6,qualityIndex = "CH")
seqPCoA(ibd_taxa,groups=groups,clusters=clusters, topTaxa=30)

## ---- message=FALSE, warning=FALSE---------------------------------------
clus.table = table(groups, clusters) 
chisq.test(clus.table) 

