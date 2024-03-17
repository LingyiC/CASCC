##' Kmeans with centers
##'
##'
##'
##' @param res The CASCC object
##' @return
##'
##' @references
##' @examples
##' library("CASCC")
##' library("cafr")
##' library("Seurat")
##' data("Data_PDAC_peng_2k") 
##' data <- Data_PDAC_peng_2k
##' # Step 0 - Step 3
##' fs.res <- CASCC.featureSelection(data, inputDataType = "well-normalized", attr.raw = NULL, exponent.max = 10, exponent.min = 2, 
##'                                  generalSeedList = NULL, mc.cores = 1, topDEGs = 10, topAttr = 50,
##'                                  removeMT.RP.ERCC = TRUE, removeNonProtein = FALSE,
##'                                  topN.DEG.as.seed = 1,
##'                                  overlapN = 10)
##' adata  <- fs.res$adata; 
##' SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
##' if (SeuratVersionCheck == 4) {
##'   data <- adata[["RNA"]]@data
##' }else if (SeuratVersionCheck == 5) {
##'   data <- adata[["RNA"]]$data
##' }
##' data <- as.matrix(data)
##' 
##' # Step 4
##' res.predictK <- predictClusterK(data, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = FALSE)
##' 
##'# results list
##'res <- list()
##'res$fs.res <- fs.res
##'res$res.predictK <- res.predictK
##'
##'# Step 5 
##'res$mainType.output <- kmeansCenterModify(res)
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

  SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
  if (SeuratVersionCheck == 4) {
    res.kmeans <- run.kmeans.fixCenters(adata[["RNA"]]@data, res$fs.res$features, attrs)
  }else if (SeuratVersionCheck == 5) {
    res.kmeans <- run.kmeans.fixCenters(adata[["RNA"]]$data, res$fs.res$features, attrs)
  }
  
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
