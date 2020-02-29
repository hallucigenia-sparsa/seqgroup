#' @title Get silhouette score for cluster labels
#'
#' @description A wrapper around vegdist that computes the silhouette score directly on these distances
#'
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param group vector with group membership
#' @param method name of dissimilarity or distance to use. See vegdist for options
#' #' @export
#'
vegdist_silhouette <- function(abundances, group, method='bray'){
  if (length(group) == length(colnames(abundances))){
    abundances <- t(abundances)
  } else if (length(group) != length(rownames(abundances))){
    stop('Group vector is not the same length as the abundance table dimensions.')
  }
  dis <- as.matrix(vegan::vegdist(abundances, method=method))
  # for every sample in a cluster, calculate silhoutte coefficient
  scores <- c()
  for (i in 1:nrow(dis)){
    row <- dis[i,]
    cluster <- group[i]
    intra_cluster_distance <- 0
    smallest_cluster_distance <- 1
    for (clus in unique(group)){
      ids <- which(group == clus)
      if (clus == cluster){
        ids <- ids[ids != i]  # distance to self should not be included
        intra_cluster_distance <- mean(row[ids])
      } else if (meandist < smallest_cluster_distance){
        smallest_cluster_distance <- mean(row[ids])
      }
    }
    if (intra_cluster_distance < smallest_cluster_distance){
      silhouette <- 1-(intra_cluster_distance / smallest_cluster_distance)

    } else if (intra_cluster_distance > smallest_cluster_distance){
      silhouette <- (smallest_cluster_distance / intra_cluster_distance) - 1
    } else {
      silhouette <- 0
    }
    scores <- c(scores, silhouette)
  }
  return(mean(scores))
}
