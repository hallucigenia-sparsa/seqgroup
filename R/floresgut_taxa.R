#' @title Gut taxon abundances from Flores et al.
#'
#' @description Data were obtained through 16S sequencing (V4 region) of faecal material.
#' @details OTU abundances were downloaded from QIITA in Nov 2020 (ID: 2151). Control and non-gut samples, samples with less than 10000 reads,
#' with a mislabel chance above 50% or no mislabel data were removed. Samples belonging to subjects with less than 6 samples
#' were also removed. Samples were rarefied to 10000 reads and OTUs were aggregated to genus level. OTUs without a known genus were summed into
#' a garbage taxon named sumUnknown. Samples were ordered first by subject and then by time (week number).
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
