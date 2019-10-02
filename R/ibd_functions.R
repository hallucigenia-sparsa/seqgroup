#' @title Functions from the PRISM cohort
#'
#' @description Function (i.e. enzyme-coding gene) abundances extracted from whole-genome shotgun data of stool samples.
#' @details Stool samples were taken from patients with Crohn's disease (CD, n=68) and ulcerative colitis (UC, n=53) as well as healthy controls (n=34).
#' The order of samples matches those of ibd_taxa. Functional assignment was carried out with the HUMAnN2 pipeline. Enzyme names as well as EC
#' numbers are given.
#' @return Function x sample matrix of size 2113 x 155
#' @docType data
#' @usage data(ibd_functions)
#' @keywords data, datasets
#' @references
#' Franzosa et al. (2019) Gut microbiome structure and metabolic activity in inflammatory bowel disease Nature Microbiology 4, 293-305
#' \href{https://www.nature.com/articles/s41564-018-0306-4}{Nature Microbiology}
#' @examples
#' data(ibd_functions)
"ibd_functions"
