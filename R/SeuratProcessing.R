##' Intial clustering
##'
##' To prepare the list of seed genes, we apply the low-complexity "Seurat"
##' workflow to the gene expression matrix.
##'
##'
##' @param data Gene expression matrix
##' @param topN By default, topN = 1. TopN ranked differentially expressed
##' genes (DEGs) of each of the clusters found in initial clusters will be used
##' as the seed gene.
##' @param resolution By default, resolution = 2. Resolution parameter in
##' \code{Seurat::FindClusters}
##' @param needNorm If the input data needs normalization. When \code{TRUE},
##' the input matrix is UMI count, we use \code{Seurat::NormalizeData} to
##' normalize.
##' @param min.diff.pct \code{Seurat::FindAllMarkers.}
##' @param topDEGs The \code{topDEGs} top-ranked DEGs of each cluster will be
##' used as features in the following analysis.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{CASCC.featureSelection}}
##' @references
##' @examples
##'
##'
##' @export


SeuratProcessing <- function(data, topN.DEG.as.seed = 1, resolution = 2, needNorm = FALSE, min.diff.pct = 0.5, topDEGs = 10){
##======================================================================
## SeuratProcessing define the seeds
##======================================================================
  # # regarding min.diff.pct:
  # if(sum(data == 0)/(ncol(data)*nrow(data)) > 0.5){ # large real dataset, the droptout rate > 0.7
  #   min.diff.pct = 0.5
  # }


  adata <- Seurat::CreateSeuratObject(counts = data)
  if (needNorm == TRUE){
    adata <- Seurat::NormalizeData(adata)
  }

  if(ncol(data) > 300){
    adata <- Seurat::FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000, verbose = F)
    adata <- Seurat::ScaleData(adata, verbose = F)
    adata <- Seurat::RunPCA(adata, verbose = F)
    adata <- Seurat::FindNeighbors(adata, verbose = F)
    adata <- Seurat::FindClusters(adata, resolution = resolution, verbose = F)
    adata <- Seurat::RunUMAP(adata, dim = 1:10, verbose = F)
  }else {
    adata <- Seurat::FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000, verbose = F)
    adata <- Seurat::ScaleData(adata, verbose = F)
    adata <- Seurat::RunPCA(adata, verbose = F)
    adata <- Seurat::RunUMAP(adata, dim = 1:10, verbose = F)
    data <- adata[["RNA"]]@data
    res_kmeans <- run.kmeans(data[Seurat::VariableFeatures(adata), ], 10)
    Seurat::Idents(adata) <- as.factor(res_kmeans$clusteringResults)
  }



## Method
  dfMarkers <-  Seurat::FindAllMarkers(adata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = min.diff.pct)
  if(dim(dfMarkers)[1] < (length(unique(Seurat::Idents(adata)))*5)){
    # Nonetheless, if the dataset lacks DEGs under this condition, for example,
    # when the average number of DEGs detected in each cluster is below 5,
    # we remove this strict restriction and include all genes in the differential expression analysis.
    # each cluster has at least ~ 5 DEGs
    dfMarkers <- Seurat::FindAllMarkers(adata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = -Inf)
    print("Warning: the dataset lacks DEGs under strict restriction")
  }
  library("dplyr")
  tmp <- dfMarkers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(n = topN.DEG.as.seed, order_by = avg_log2FC)
  markers <- tmp$gene
  markers <- unique(markers)
  # markers <- intersect(markers, rownames(data))




  tmp <- dfMarkers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(n = topDEGs, order_by = avg_log2FC)
  DEGs <- tmp$gene
  DEGs <- unique(DEGs)
  # DEGs <- intersect(DEGs, rownames(data))


  res <- list()
  res$adata <- adata
  res$markers <- markers
  res$DEGs <- DEGs
  res$dfMarkers <- dfMarkers
  return(res)
}

