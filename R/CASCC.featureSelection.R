##' CASCC Feature Selection
##'
##' This function help select features for clustering.
##'
##'
##' @param data Gene expression matrix.
##' @param inputDataType Input data type. Default is "well-normalized".
##' Options: "well-normalized", "needLog2" or "UMI". If the input matrix is
##' normalized and log-transformed, use "well-normalized". If the input matrix
##' is normalized but not log-transformed, use "needLog2". If the input matrix
##' is UMI-count, use "UMI".
##' @param attr.raw Raw attractor object. Default is NULL.  If users already
##' obtained attractor scanning results \code{CASCC.output$fs.res$attrs.raw},
##' they can provide the list of attractors to skip the scanning step.
##' @param exponent.max Maximum exponent for adaptive attractor. Default is 10.
##' @param exponent.min Minimum exponent for adaptive attractor. Default is 2.
##' @param generalSeedList Seed genes set by users. Default is NULL.
##' @param mc.cores Number of cores to use.
##' @param topDEGs Number of top DEGs in each cluster selected as features.
##' Default is 10.
##' @param topAttr Number of top genes in each attractor selected as features.
##' Default is 50.
##' @param removeMT.RP.ERCC If remove mitochondrial genes, ribosomal genes and
##' spike-in genes from the analysis. Default is TRUE.
##' @param topN.DEG.as.seed The Number of top DEGs to use as seeds. Default is
##' 1.
##' @param removeNonProtein If remove non-protein coding genes.  Default is
##' FALSE.
##' @param overlapN Two attractors are considered duplicated if their top 50 genes
##' overlap by N genes. The default is 10.
##' @return \item{features}{Selected features.} \item{adata}{Seurat object
##' after gene filtering.} \item{dfMarkers}{Differential expression genes for
##' initial clusters.} \item{markers}{Seeds for attractors.}
##' \item{attrs.raw}{Raw results from \code{link{scanSeeds.adaptive}.}}
##' \item{attrs}{A list of attractors in which identical attractors were
##' removed.} \item{finalAttrs}{A list of attractors in which duplicated and
##' identical attractors were removed.}
##' @author Lingyi Cai
##' @seealso \code{\link{filterGenes}}
##' @references
##' @examples
##'
##'
##' @export

CASCC.featureSelection <- function(data, inputDataType = "well-normalized", attr.raw = NULL, exponent.max = 10, exponent.min = 2, generalSeedList = NULL, mc.cores = NULL, topDEGs = 10, topAttr = 50,
                                   removeMT.RP.ERCC = TRUE, topN.DEG.as.seed = 1, removeNonProtein = FALSE, overlapN = 10,
                                   attractorModify.percentage = 0.01){
  set.seed(1)
  resList <- list()

  # preprocess-Filter genes
  data <- filterGenes(data, min.cells = 0, removeMT.RP.ERCC = removeMT.RP.ERCC, removeNonProtein = removeNonProtein)

  if(inputDataType == "well-normalized"){
    needNorm = FALSE
  }else if(inputDataType == "needLog2"){
    needNorm = FALSE
    data = log2(data + 1)
  }else if(inputDataType == "UMI"){
    needNorm = TRUE
  }

  # Step 1: initial clustering
  seurat.res <- SeuratProcessing(data, needNorm = needNorm, topDEGs = topDEGs, topN.DEG.as.seed = topN.DEG.as.seed)
  markers <- seurat.res$markers
  DEGs <- seurat.res$DEGs
  adata <- seurat.res$adata
  data <- adata[["RNA"]]@data

  if(is.null(generalSeedList) == FALSE){
    if(generalSeedList == "generalSeedList"){
      generalSeedList = c("AMBP", #ductal cell 1
                          "KRT19", #ductal cell 2
                          "PRSS1", #acinar
                          "CHGB", # endocrine
                          "RGS5", # stellate
                          "LUM", # fibroblast
                          "CDH5", # endothelial
                          "AIF1",# macrophage
                          "CD3D", # T cell
                          "MS4A1") # B cell
      generalSeedList <- intersect(rownames(data), generalSeedList)
      markers <- unique(c(generalSeedList, markers))
    }else{
      markers <- unique(c(generalSeedList, markers))
      markers <- intersect(rownames(data), markers)
    }
  }

  # Step 2: adaptive attractor
  if(is.null(attr.raw) == TRUE){
    res <- scanSeeds.adaptive(data, markers, exponent.max = exponent.max, exponent.min = exponent.min, mc.cores = mc.cores)
  }else {
    res <- attr.raw
  }


  # Step 3: selected features
  # attrs <- collapseAttractorList(res) ## 50 50
  attrs <- collapseAttractorList(res, NumTopFeature=50) # remove the identical attractors 
  features <- lapply(attrs, function(x){names(x[1:topAttr])})
  # features <- lapply(attrs, function(x){names(x[1: min(round(nrow(data)*0.005), 50)])})
  features <- unique(unlist(features))
  features <- unique(c(features, DEGs))


  resList$features <- features #topN attrs + DEGs
  resList$adata <- adata
  resList$dfMarkers <- seurat.res$dfMarkers
  resList$markers <- markers #seeds
  resList$attrs.raw <- res
  resList$attrs <- attrs # only remove identical attractors


  # Generate final Attractors

  # ./R/collapseAttractorList.R
  # remove similar attractors
  # remove the duplicated attractors 
  attrs <- collapseAttractorList(res, NumTopFeature=50, overlapN = overlapN) # remove the duplicated attractors 
  ## backup 
  attrs <-  attractorModify(attrs, adata, expressionLevel = 1, percentage = attractorModify.percentage) # default attractorModify.percentage is 0, equals we skip this step

  resList$finalAttrs <- attrs # removed duplciated attractors

  return(resList)
}
