#' @title Metadata for the PRISM cohort
#'
#' @description Metadata include SRA identifier, age, diagnosis, fecal calprotectin, antibiotic and immunosuppressant usage as well as treatment with mesalamine or steroids.
#' @details Stool samples were taken from patients with Crohn's disease (CD, n=68) and ulcerative colitis (UC, n=53) as well as healthy controls (n=34).
#' CD/UC vs healthy status is stored in 'Diagnosis'. The order of samples matches those of ibd_taxa.
#' @return Sample x metadata dataframe of size 155 x 8
#' @docType data
#' @usage data(ibd_metadata)
#' @keywords data, datasets
#' @references
#' Franzosa et al. (2019) Gut microbiome structure and metabolic activity in inflammatory bowel disease Nature Microbiology 4, 293-305
#' \href{https://www.nature.com/articles/s41564-018-0306-4}{Nature Microbiology}
#' @examples
#' data(ibd_metadata)
#' # convert to numeric
#' ibd_metadata$Fecal.Calprotectin=as.numeric(ibd_metadata$Fecal.Calprotectin)
#' ibd_metadata$Age=as.numeric(ibd_metadata$Age)
#' summary(ibd_metadata)
"ibd_metadata"
