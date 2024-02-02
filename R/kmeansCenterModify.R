##' Kmeans with centers
##'
##'
##'
##' @param res The CASCC project
##' @return
##'
##' @references
##' @examples
##'
##'
##' @export

kmeansCenterModify <- function(res){
  output <- list()
  adata <- res$res.predictK$adata


  # ## scenario 1
  # attrs <- collapseAttractorList(res$fs.res$attrs.raw, NumTopFeature=50, overlapN = 50, verbose = F) #remove identical attractors
  # # ./R/collapseAttractorList.R
  # # remove duplicate attractors
  # ## scenario 2
  # attrs <- collapseAttractorList(res$fs.res$attrs.raw, NumTopFeature=50, overlapN = 10)
  # ## scenario 3
  # attrs <-  attractorModify(attrs, res$fs.res$adata, expressionLevel = 1, percentage = 0.01)

  attrs <- res$fs.res$finalAttrs
  # attrs <- res$fs.res$attrs
  attrs <- attrs[which(rank(unlist(lapply(attrs, function(x){x[1] - x[2]}))) <= res$res.predictK$clusterK)]

  # attrs <- attrs[which(rank(unlist(lapply(attrs, function(x){-x[5]}))) <= max(3, res$res.predictK$clusterK))]

  res.kmeans <- run.kmeans.fixCenters(adata[["RNA"]]@data, res$fs.res$features, attrs)
  clusteringResults <- res.kmeans$clusteringResults


  output$clusteringResults <- clusteringResults
  output$attrs.center <- attrs
  return(output)
}

##' Kmeans with centers Using attractors as centers
##'
##'
##' @param data.full The full data matrix
##' @param features The features
##' @param attrs The attractors
##' @return
##'
##' @references
##' @examples
##'
##'
##' @export

### Using attractors as centers
run.kmeans.fixCenters <- function(data.full, features, attrs){
  set.seed(1)
  res <- list()
  # attrs <- collapseAttractorList(attrs.raw, NumTopFeature=50, overlapN = 10)
  data.full <- as.matrix(data.full)
  centers <- unique(names(unlist(lapply(attrs, function(x){which.max(colMeans(data.full[names(x[1:5]), ]))}))))
  res$obj <- stats::kmeans(t(data.full[features, ]), t(data.full[features, centers]), iter.max = 1000, nstart = 100)
  res$clusteringResults <- as.factor(res$obj$cluster)

  ## remove those attractors do not have many cells
  # indices <- as.numeric(which(table(res$clusteringResults) < ncol(data.full)*0.001))

  # if ( length(indices) != 0) {
  #   attrs <- attrs[-indices]
  #   centers <- unique(names(unlist(lapply(attrs, function(x){which.max(colMeans(data.full[names(x[1:5]), ]))}))))
  #
  #   res$obj <- stats::kmeans(t(data.full[features, ]), t(data.full[features, centers]), iter.max = 1000, nstart = 100)
  #   res$clusteringResults <- as.factor(res$obj$cluster)
  #   res$centers <- centers
  # }

  return(res)
}

##=== Kmeans
# > initial.centers <- names(unlist(lapply(attrs, function(x){which.max(apply(data[names(x[1:50]), ], 2, mean))})))
##' Run kmeans
##'
##'
##' @param data.lognormcount The log normalized count matrix
##' @param centers The number of centers
##' @return
##'
##' @references
##' @examples
##'
##'
##' @export
run.kmeans <- function(data.lognormcount, centers){
  set.seed(1)
  res <- list()
  data.lognormcount <- as.matrix(data.lognormcount)
  res$obj <- stats::kmeans(t(data.lognormcount), centers, iter.max = 1000, nstart = 100)
  res$clusteringResults <- as.factor(res$obj$cluster)
  return(res)
}
