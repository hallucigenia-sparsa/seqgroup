#' @title Significance of network intersection
#'
#' @description Builds groups of networks randomly and counts their number of intersection edges (intersection size).
#' Each random network is built by randomly sampling from the total number of possible edges.
#' The p-value is computed parameter-free by counting the number of times the intersection size of
#' random networks is smaller than the intersection size in the observed networks and dividing it
#' by the number of iterations. Please see \href{https://github.com/ramellose/anuran}{https://github.com/ramellose/anuran}
#' for a dedicated tool that assesses the significance of network intersections.
#'
#' @param rep.num number of inferred networks
#' @param taxon.num number of taxa in the matrix
#' @param avg.network.size average number of edges in inferred networks
#' @param inter.size number of edges in observed intersection network
#' @param iter number of iterations for random network group construction
#' @param directed if directed, the possible edge number is taxon.num * (taxon.num-1), else it is (taxon.num * (taxon.num-1))/2
#' @param loops if directed, count self-arcs (so compute possible edge number as taxon.num * taxon.num)
#' @param hist.rand if TRUE, plot the histogram of random intersection sizes
#' @return p-value
#' @examples
#' networkIntersectSig(rep.num=3,taxon.num=8,avg.network.size=12,inter.size=10,directed=TRUE)
#' @export
networkIntersectSig<-function(rep.num=5, taxon.num=20, avg.network.size=0,inter.size=0, iter=1000, directed=FALSE, loops=FALSE, hist.rand=FALSE){
  rand.intersections=c()
  # total possible edge number for given number of taxa
  total.edge.num=NA
  if(directed){
    if(loops){
      total.edge.num=taxon.num*taxon.num
    }else{
      total.edge.num=taxon.num*(taxon.num-1)
    }
  }else{
    total.edge.num=(taxon.num*(taxon.num-1))/2
  }
  bigger.than.obs=0
  # repeat the simulation of the random network
  for(i in 1:iter){
    intersection=c()
    for(rep in 1:rep.num){
      # assume each edge is equally likely, take a random sub-set of the average size of the inferred networks
      edges.rand=sample(1:total.edge.num)[1:avg.network.size]
      if(rep>1){
        intersection=intersect(intersection,edges.rand)
      }else{
        intersection=edges.rand
      }
    }
    rand.intersections=c(rand.intersections,length(intersection))
    # record the number of times the random intersection network is bigger than the observed intersection network
    if(length(intersection)>=inter.size){
      bigger.than.obs=bigger.than.obs+1
    }
  }
  print(paste("Mean random intersections: ",mean(rand.intersections)))
  print(paste("Standard deviation random intersections: ",sd(rand.intersections)))
  if(hist.rand){
    hist(rand.intersections,main="Histogram of random intersection sizes", xlab="Random intersection size",xlim=c(0,inter.size))
    abline(v=inter.size,col="red")
  }
  pval=(bigger.than.obs+1)/(iter+1)
  return(pval)
}
