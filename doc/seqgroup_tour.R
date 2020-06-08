## ----setup, include=FALSE-----------------------------------------------------
# Global options
library(knitr)
opts_chunk$set(fig.path="figure_seqgroup_tour/")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(seqgroup)

## ---- message=FALSE, warning=FALSE--------------------------------------------
data("ibd_taxa") # taxon x sample matrix
data("ibd_lineages") # lineage x taxon matrix
data("ibd_metadata") # metadata x sample data frame
data("ibd_functions") # functions x sample matrix
summary(ibd_metadata)

## ---- message=FALSE, warning=FALSE--------------------------------------------
groups=as.vector(ibd_metadata$Diagnosis)
ibd.groups=groups
ibd.groups[ibd.groups=="UC"]="IBD"
ibd.groups[ibd.groups=="CD"]="IBD"

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
compareGroups(ibd_taxa,groups=groups,property="alpha",pvalViz = TRUE)
compareGroups(ibd_taxa,groups=groups,property="beta",pvalViz = TRUE)
compareGroups(ibd_functions,groups=groups,property="beta",pvalViz = TRUE)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 3, fig.width = 3----
compareGroups(ibd_taxa,groups=groups,property="beta",method="DM")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
control.indices=which(groups=="Control")
groupBarplot(ibd_taxa[,control.indices],topTaxa = 15,legend.shift = 1.5, main="Control")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
groupBarplot(ibd_taxa,topTaxa = 15,groups=groups, randSampleNum = 10, legend.shift = 1.5, main="Taxon composition")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
ibd.indices=which(ibd.groups=="IBD")
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
ibd.groups.nafree=ibd.groups[indices.to.keep]
taxon="Faecalibacterium_prausnitzii"
compareDistribs(ibd_taxa.nafree,taxon=taxon,groups=ibd.groups.nafree,group1="IBD",group2="Control")

## ---- message=FALSE, results='hide', warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
clusters=findClusters(ibd_taxa,k=NA,maxK=5,qualityIndex = "CH")
seqPCoA(ibd_taxa,groups=groups,clusters=clusters, topTaxa=30)

## ---- message=FALSE, warning=FALSE--------------------------------------------
clus.table = table(ibd.groups, clusters) 
chisq.test(clus.table) 

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 10, fig.width = 10----
def.par=par(no.readonly = TRUE) # save previous par settings
layout(matrix(c(1,3,2,4),2,2,byrow = TRUE), c(2,3,2,3), TRUE) # define space for 4 figures
ibd.network=buildNetwork(ibd_taxa[indices.prev,ibd.indices], lineages=ibd_lineages[indices.prev,], method=c("spearman","bray"), repNum=0, nameLevel="species") # run CoNet on IBD samples
control.network=buildNetwork(ibd_taxa[indices.prev,control.indices], lineages=ibd_lineages[indices.prev,], method=c("spearman","bray"), repNum=0, nameLevel="species") # run CoNet on control samples
plot(ibd.network,main="IBD network")
plot(control.network, main="Control network")
par(def.par)

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
control.aracne=buildNetwork(ibd_taxa[indices.prev,control.indices],lineages=ibd_lineages[indices.prev,],method="aracne",repNum=100)
plot(control.aracne, main="Control network inferred with ARACNE")
# The node label is assigned on genus level here; unlabeled nodes do not have a genus assigned.

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 6, fig.width = 6----
# we binarise the group assignment to append it to the abundances
ibd.status=ibd.groups
ibd.status[ibd.groups=="Control"]=0
ibd.status[ibd.groups=="IBD"]=1
ibd_taxa_extended=rbind(ibd_taxa[indices.prev,],as.numeric(ibd.status))
rownames(ibd_taxa_extended)=c(rownames(ibd_taxa)[indices.prev],"ibd.status") # append row name for IBD status
ibd_lineages_extended=rbind(ibd_lineages[indices.prev,],rep("IBD_status",7)) # entry in lineage matrix needed to label node
rownames(ibd_lineages_extended)=c(rownames(ibd_lineages)[indices.prev],"ibd.status") # append row name for IBD status
network=buildNetwork(ibd_taxa_extended, lineages=ibd_lineages_extended, method=c("spearman","bray"), repNum=0, nameLevel="species") # run CoNet on all samples
plot(network,main="Global network")

## ---- message=FALSE, warning=FALSE, fig.align = "center", fig.height = 10, fig.width = 10----
prev=150 # minimum occurrence for functions
ibd_functions_incidence=ibd_functions
ibd_functions_incidence[ibd_functions_incidence>0]=1
rowSums=rowSums(ibd_functions_incidence)
indices.functions.prev=which(rowSums>=prev)
# CoNet expects metadata to have samples as rows, therefore functions have to be transposed
bip.network=barebonesCoNet(abundances=ibd_taxa[indices.prev,],metadata=t(ibd_functions[indices.functions.prev,]),min.occ=0,methods="pearson",T.up=0.7,T.down=-0.7) 
plot(bip.network,main="Global species-function network")

