#' @title Gut taxon abundances from Flores et al.
#'
#' @description Data were obtained through 16S sequencing of the V4 region.
#' @details OTU abundances were downloaded from QIITA in Nov 2020 (ID: 2151). Control and non-gut samples, samples with less than 10000 reads,
#' mislabel chance above 50% or no mislabel data were removed. Samples belonging to subjects with less than 6 samples taken
#' were also removed. Samples were rarefied to 10000 reads and OTUs were aggregated to genus level. OTUs without a known genus were summed into
#' a garbage taxon named sumUnknown. The resulting genus count matrix has 323 taxa and 637 samples.
#' @return Taxon x sample matrix of size 323 x 637
#' @docType data
#' @usage data(floresgut_taxa)
#' @keywords data, datasets
#' @references
#' Flores et al. (2014) Temporal variability is a personalized feature of the human microbiome Genome Biology 15, 531
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0531-y}{Genome Biology}
#' @examples
#' data(floresgut_taxa)
"floresgut_taxa"
