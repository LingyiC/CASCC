##' Predict the number of clusters
##'
##'
##'
##' @param data A matrix of gene expression data. Rows are genes and columns
##' are cells.
##' @param fs.res A list of feature selection results generated by
##' \code{CASCC.featureSelection} function.
##' @param method The cluster analysis method to be used. Default is "ward.D2".
##' @param index The index to be calculated. Default is "silhouette".
##' @param min.nc_fix Whether to fix the minimum number of clusters. Default is
##' FALSE.
##' @param if_PCA_input Whether to use PCA to predict the number of clusters. Default is
##' FALSE. Then use UMAP 2D embedding
##' @param nPCA The number of PCA. Default is NULL
##' @return \item{clusterK}{The estimated number of clusters.}
##' \item{adata}{Seurat object after CASCC feature selection.}
##' @author Lingyi Cai
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
##' @export


predictClusterK <- function(data, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = TRUE, if_PCA_input = FALSE, nPCA = NULL){
  res.predictK <- list()
  features <- fs.res$features

  adata <- Seurat::CreateSeuratObject(counts = data) # this step is the full

  SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
  if (SeuratVersionCheck == 4) {
    adata@assays$RNA@var.features <- features
  }else if (SeuratVersionCheck == 5) {
    adata[["RNA"]]$data <- data
    VariableFeatures(adata) <- features
  }
  adata <- Seurat::ScaleData(adata, verbose = F)
  adata <- Seurat::RunPCA(adata, verbose = F)
  adata <- Seurat:: RunUMAP(adata, dim = 1:10, verbose = F)
  # source("/drive_T7/2021Research/MicroE/TASK/2023/0125_CASCC_single_function/collapseAttractorList.R")

  finalAttrs <- fs.res$finalAttrs
  # We modify the list of attractors
  # according to two criteria: First, we filter out any attractor in which the score difference between
  # the 1 st and 2 nd top genes is greater than 0.2, indicating that the attractor is dominated by the top
  # gene. Second, we keep only attractors where the score of the 5 th ranked gene is at least 0.4,
  # which suggests a sufficiently strong co-expression mechanism involving at least five genes.
  selected <- unlist(lapply(finalAttrs, function(x) {x[5] >= 0.4 && (x[1]-x[2]) < 0.2}) )
  attrs.veryStable <- finalAttrs[selected]

  # max.nc equals the number of attractors after removing the duplicated attratcors
  if(if_PCA_input == FALSE){
    if(min.nc_fix == TRUE){
      res.nbclust <- NbClust::NbClust(adata@reductions$umap@cell.embeddings,
                                      min.nc = 3, max.nc = max(3, length(fs.res$finalAttrs)), #fs.res$attrs
                                      method = method, index =index)
    }else {
      # default options
      res.nbclust <- NbClust::NbClust(adata@reductions$umap@cell.embeddings,
                                      min.nc = max(3, length(attrs.veryStable)), max.nc = max(3, length(fs.res$finalAttrs)), #attr.res$attractorSignatureList
                                      method = method, index =index)
    }

  }else if(if_PCA_input == TRUE){
    if(min.nc_fix == TRUE){
      res.nbclust <- NbClust::NbClust(adata@reductions$pca@cell.embeddings[, 1:nPCA],
                                      min.nc = 3, max.nc = max(3, length(fs.res$finalAttrs)), #fs.res$attrs
                                      method = method, index =index)
    }else {
      # default options
      res.nbclust <- NbClust::NbClust(adata@reductions$pca@cell.embeddings[, 1:nPCA],
                                      min.nc = max(3, length(attrs.veryStable)), max.nc = max(3, length(fs.res$finalAttrs)), #attr.res$attractorSignatureList
                                      method = method, index =index)
    }
  }

  res.predictK$clusterK <- res.nbclust$Best.nc["Number_clusters"]
  res.predictK$umap.cell.embeddings <- adata@reductions$umap@cell.embeddings
  res.predictK$adata <- adata
  return(res.predictK)
}



# predictClusterK_2 <- function(data, attr.res, features){
#   res.predictK <- list()
#   adata <- CreateSeuratObject(counts = data, verbose = F)
#   adata@assays$RNA@var.features <- features
#   adata <- ScaleData(adata, verbose = F)
#   adata <- RunPCA(adata, verbose = F)
#   adata <- RunUMAP(adata, dim = 1:10, verbose = F)
#   source("/drive_T7/2021Research/MicroE/TASK/2023/0125_CASCC_single_function/collapseAttractorList.R")
#   attrs <- collapseAttractorList(res, NumTopFeature=50, overlapN = 5)
#   res.nbclust <- NbClust::NbClust(adata@reductions$umap@cell.embeddings,
#                   min.nc = max(2, length(attrs)), max.nc = length(collapseAttractorList(res, NumTopFeature=50, overlapN = 50)),
#                   method = "ward.D2", index =c("ch"))

#   res.predictK$umap.cell.embeddings <- adata@reductions$umap@cell.embeddings
#   res.predictK$clusterK <- res.nbclust$Best.nc["Number_clusters"]
#   return(res.predictK)
# }
