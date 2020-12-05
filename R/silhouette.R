#' @title Get silhouette score for a clustering result
#'
#' @description The silhouette score is a cluster qualiy index that assessses both cluster cohesion
#' (how close points within the same cluster are to each other) and cluster separation
#' (how well different clusters are separated) for a given clustering result (encoded in the group membership vector).
#' The silhouette score is bounded by -1 (worst value) and 1 (best value). The silhouette score s of a point i is defined as:
#' s(i)=b(i)-a(i)/max(b(i),max(i)), where a(i) is the mean distance of this point to all other points within the same cluster
#' and b(i) is the smallest mean distance to data points in the nearest neighboring cluster.
#' The silhouette score of the entire clustering is computed as the mean of silhouette scores across all points.
#' This function wraps vegdist to compute distances/dissimilarities.
#'
#' @param abundances a matrix with taxa as rows and samples as columns
#' @param groups vector with group memberships
#' @param method name of dissimilarity or distance to use. See vegdist for options
#' @export
#'
silhouette <- function(abundances, groups, method='bray'){
  if (length(groups) == length(colnames(abundances))){
    abundances <- t(abundances)
  } else if (length(groups) != length(rownames(abundances))){
    stop('Group vector is not the same length as the abundance table dimensions.')
  }
  #print(dim(abundances))
  dis <- as.matrix(vegdist(abundances, method=method))
  # for every sample in a cluster, calculate silhoutte coefficient
  scores <- c()
  empty=FALSE
  for (i in 1:nrow(dis)){
    row <- dis[i,]
    cluster <- groups[i]
    intra_cluster_distance <- 0
    smallest_cluster_distance <- 1
    for (clus in unique(groups)){
      ids <- which(groups == clus)
      if (clus == cluster){
        ids <- ids[ids != i]  # distance to self should not be included
        # Manual merge of treatment of zero-member clusters with alternative
        # <<<<<<< HEAD
        #        if(length(ids)==0){
        #          warning(paste("Cluster",cluster," of",length(unique(groups)),"clusters is empty! Silhouette is set to NA"))
        #          empty=TRUE
        #        }else{
        #          #print(paste("ids=",ids))
        #          intra_cluster_distance <- mean(row[ids])
        #        }
        #        # fix: replaced meandist below by intra_cluster_distance
        #      } else if (intra_cluster_distance < smallest_cluster_distance){
        #        smallest_cluster_distance <- mean(row[ids])
        # =======
        # only done if there are members in the cluster: intra_cluster_distance <- mean(row[ids])
        if(length(ids)==0){
          warning(paste("Cluster",cluster," of",length(unique(groups)),"clusters is empty! Silhouette is set to NA"))
          empty=TRUE
        }else{
          #print(paste("ids=",ids))
          intra_cluster_distance <- mean(row[ids])
        }
      } else {
        inter_cluster_distance <- mean(row[ids])
        if (inter_cluster_distance < smallest_cluster_distance){
          smallest_cluster_distance <- inter_cluster_distance
        }
# >>>>>>> 0d6c5931d0584caae022c3a1319236ef0c91ed4c
      }
    } # loop clusters
    if (intra_cluster_distance < smallest_cluster_distance){
      silhouette <- 1-(intra_cluster_distance / smallest_cluster_distance)

    } else if (intra_cluster_distance > smallest_cluster_distance){
      silhouette <- (smallest_cluster_distance / intra_cluster_distance) - 1
    } else {
      silhouette <- 0
    }
    if(empty){
      silhouette=NA
    }
    #print(paste("silhouette",silhouette))
    scores <- c(scores, silhouette)
  } # loop rows of dissimilarity matrix
  #print(scores)
  return(mean(scores))
}
