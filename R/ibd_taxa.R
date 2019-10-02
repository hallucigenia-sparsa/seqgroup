#' @title IBD and control taxon abundances from the PRISM cohort
#'
#' @description Species abundances extracted from whole-genome shotgun data of stool samples.
#' @details Samples include patients with Crohn's disease (CD, n=68) and ulcerative colitis (UC, n=53) as well as healthy controls (n=34).
#' Species relative abundances were computed with MetaPhlAn2 and species below or equal to 0.1% abundance in at least five samples were excluded.
#' @return Species x sample matrix of size 201 x 155
#' @docType data
#' @usage data(ibd_taxa)
#' @keywords data, datasets
#' @references
#' Franzosa et al. (2019) Gut microbiome structure and metabolic activity in inflammatory bowel disease Nature Microbiology 4, 293-305
#' \href{https://www.nature.com/articles/s41564-018-0306-4}{Nature Microbiology}
#' @examples
#' data(ibd_taxa)
"ibd_taxa"
