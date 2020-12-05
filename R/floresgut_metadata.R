#' @title Gut sample metadata from Flores et al.
#'
#' @description Metadata are selected from the sample information file in QIITA (ID: 2151) and include age, bmi, diet, gender,
#' host_subject_id, ibd, menstruation_disturbance, raceethnicity, sickness_disturbance, smokecigarettes and weeknumber.
#' Metadata is a data frame with samples as rows and metadata items as columns. Samples are in the same order as the samples in
#' floresgut_taxa.
#' @details Note that some metadata items contain missing values.
#' @return Sample x metadata data frame of size 637 x 11
#' @docType data
#' @usage data(floresgut_metadata)
#' @keywords data, datasets
#' @references
#' Flores et al. (2014) Temporal variability is a personalized feature of the human microbiome Genome Biology 15, 531
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0531-y}{Genome Biology}
#' @examples
#' data(floresgut_metadata)
"floresgut_metadata"
