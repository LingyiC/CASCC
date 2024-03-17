##' Attractor Plot
##'
##' This function colors singel cell on UMAP according to the expression levels
##' of co-expressed signatures.
##'
##'
##' @param seuratObj A Seurat object. For example,
##' \code{adata <- CASCC.output$res.predictK$adata}; \code{Idents(adata) <- CASCC.output$mainType.output$clusteringResults}. 
##' @param attr.list A list of attractors to plot. For example, 
##' \code{CASCC.output$mainType.output$attrs.center},
##' \code{CASCC.output$fs.res$finalAttrs}.
##' @param topN Specify the topN of each attractor to plot.
##' @param n.col Specify the number of columns in which the plots should be
##' arranged when they are displayed together.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{attractorDimplot}}
##' @references
##' @examples
##' library(CASCC)
##' library(Seurat)
##'
##' # load the toy dataset
##' data("Res_PDAC_peng_2k")
##' res <- Res_PDAC_peng_2k
##'
##' # assign clustering results to adata
##' adata <- res$res.predictK$adata # feature-selected UMAP
##' Idents(adata) <- res$mainType.output$clusteringResults # this step is necessary 
##' 
##' # plot the attractor plot
##' attractorPlot(adata, res$mainType.output$attrs.center, n.col = 4)
##'
##' @export

## a list of attractors
attractorPlot <- function(seuratObj, attr.list, topN = 5, n.col = 4){
  SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
  library("tidyverse")
  library("RColorBrewer")
  library("Seurat")
  figures <- list()
  for(i in seq(1, length(attr.list))){
    attractor <- names(attr.list[[i]])[1:topN]

  if (SeuratVersionCheck == 4) {
    mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]@data[attractor, ]), na.rm = TRUE)
  }else if (SeuratVersionCheck == 5) {
    mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]$data[attractor, ]), na.rm = TRUE)
  }

    attr_name <- paste0(names(attr.list[[i]])[1], ("_attractor"))
    seuratObj@meta.data[, attr_name] <- mean.exp
    figures[[i]] <- FeaturePlot(object = seuratObj, features = attr_name)+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  }

  library("gridExtra")
  figure <- do.call("grid.arrange", c(figures, ncol=n.col))
  return(figure)
}


##' Attractor Dimplot
##'
##' This function plot the CASCC UMAP reuslts. Each dot will be labeled with
##' the attractor it belongs to.
##'
##'
##' @param seuratObj A Seurat object. For example,
##' \code{CASCC.output$res.predictK$adata}; \code{Idents(adata) <- CASCC.output$mainType.output$clusteringResults}. 
##' @param attr.list A list of attractors to plot. For example,
##' \code{CASCC.output$mainType.output$attrs.center}.
##' @param label Whether to label attractor on the plot.
##' @param cols A vector of colors to use for each attractor. If not provided, use the default color.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{attractorPlot}}
##' @references
##' @examples
##' library(CASCC)
##' library(Seurat)
##'
##' # load the toy dataset
##' data("Res_PDAC_peng_2k")
##' res <- Res_PDAC_peng_2k
##'
##' # assign clustering results to adata
##' adata <- res$res.predictK$adata # feature-selected UMAP
##' Idents(adata) <- res$mainType.output$clusteringResults # this step is necessary 
##' 
##' # plot the attractor dimplot
##' attractorDimplot(adata, res$mainType.output$attrs.center)
##' @export

attractorDimplot <- function(seuratObj, attr.list, label = F, cols = NULL){

  SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
  ## Based on the kmeans clustering results, we assign cluster to an attractor.
  attr.list <- attr.list[which(rank(unlist(lapply(attr.list, function(x){x[1] - x[2]}))) <= length(unique(Idents(seuratObj))))]
  library("tidyverse")
  library("RColorBrewer")
  library("Seurat")
  library("ggplot2")
  attr.meta <- matrix(data = NA, nrow = length(attr.list), ncol = ncol(seuratObj))
  attr_names <- c()
  for(i in seq(1, length(attr.list))){
    attractor <- names(attr.list[[i]])[1:5] # expression of top 5 genes of each attractor

    if (SeuratVersionCheck == 4) {
      mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]@data[attractor, ]), na.rm = TRUE)
    }else if (SeuratVersionCheck == 5) {
      mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]$data[attractor, ]), na.rm = TRUE)
    }

    attr_name <- paste0(names(attr.list[[i]])[1], ("_attractor"))
    attr_names <- c(attr_names, attr_name)
    attr.meta[i, ] <- mean.exp
  }

  colnames(attr.meta) <- colnames(seuratObj)
  rownames(attr.meta) <- attr_names

  clusteringResults <- Idents(seuratObj)
  attractor_tags <- c()
  ## find based on highest expression
  # for (i in levels(clusteringResults)) {
  #   if(length(which(clusteringResults == i)) > 1){
  #     tmp1 <- rowMeans(attr.meta[, names(clusteringResults[which(clusteringResults == i)])]) # cluster_i metagene
  #   }else{
  #     tmp1 <- attr.meta[, names(clusteringResults[which(clusteringResults == i)])]
  #   }
  #
  #   tmp2 <- rowMeans(attr.meta[, names(clusteringResults[which(clusteringResults %in% setdiff(levels(clusteringResults), i))])])  # cluster_others metagene
  #   # attractor_tags <- c(attractor_tags, names(which.max(tmp1 - tmp2)))
  #   attractor_tags <- c(attractor_tags, names(which.max((tmp1 - tmp2)/(tmp2 + 0.0001))))
  # }
  ## find based on center order
  # print(identical(attractor_tags, attr_names))
  attractor_tags <- attr_names
  levels(clusteringResults) <- attractor_tags
  Idents(seuratObj) <- clusteringResults
  figure <- DimPlot(seuratObj, label = label, cols = cols)
  return(figure)
}


