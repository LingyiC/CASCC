##' Attractor Plot
##' 
##' This function colors singel cell on UMAP according to the expression levels
##' of co-expressed signatures.
##' 
##' 
##' @param seuratObj A Seurat object. For example,
##' \code{CASCC.output$res.predictK$adata}.
##' @param attr.list A list of attractors to plot. For example,
##' \code{CASCC.output$mainType.output$attrs.center},
##' \code{CASCC.output$fs.res$finalAttrs}.
##' @param n.col Specify the number of columns in which the plots should be
##' arranged when they are displayed together.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{attractorDimplot}}
##' @references
##' @examples
##' 
##' 
##' @export

## a list of attractors
attractorPlot <- function(seuratObj, attr.list, n.col = 4){
  library("tidyverse")
  library("RColorBrewer")
  library("Seurat")
  figures <- list()
  for(i in seq(1, length(attr.list))){
    attractor <- names(attr.list[[i]])[1:5]
    mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]@data[attractor, ]), na.rm = TRUE)
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
##' \code{CASCC.output$res.predictK$adata}.
##' @param attr.list A list of attractors to plot. For example,
##' \code{CASCC.output$mainType.output$attrs.center}.
##' @param label Whether to label attractor on the plot.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{attractorPlot}}
##' @references
##' @examples
##' 
##' 
##' @export

attractorDimplot <- function(seuratObj, attr.list, label = F, cols = NULL){
  ## Based on the kmeans clustering results, we assign cluster to an attractor.
  library("tidyverse")
  library("RColorBrewer")
  library("Seurat")
  library("ggplot2")
  attr.meta <- matrix(data = NA, nrow = length(attr.list), ncol = ncol(seuratObj))
  attr_names <- c()
  for(i in seq(1, length(attr.list))){
    attractor <- names(attr.list[[i]])[1:5]
    mean.exp <- colMeans(x = as.matrix(seuratObj[["RNA"]]@data[attractor, ]), na.rm = TRUE)
    attr_name <- paste0(names(attr.list[[i]])[1], ("_attractor"))
    attr_names <- c(attr_names, attr_name)
    attr.meta[i, ] <- mean.exp
  }

  colnames(attr.meta) <- colnames(seuratObj)
  rownames(attr.meta) <- attr_names

  clusteringResults <- Idents(seuratObj)
  attractor_tags <- c()
  for (i in levels(clusteringResults)) {
    if(length(which(clusteringResults == i)) > 1){
      tmp1 <- rowMeans(attr.meta[, names(clusteringResults[which(clusteringResults == i)])]) # cluster_i metagene
    }else{
      tmp1 <- attr.meta[, names(clusteringResults[which(clusteringResults == i)])]
    }

    tmp2 <- rowMeans(attr.meta[, names(clusteringResults[which(clusteringResults %in% setdiff(levels(clusteringResults), i))])])  # cluster_others metagene
    # attractor_tags <- c(attractor_tags, names(which.max(tmp1 - tmp2)))
    attractor_tags <- c(attractor_tags, names(which.max((tmp1 - tmp2)/(tmp2 + 0.0001))))
  }

  levels(clusteringResults) <- attractor_tags
  Idents(seuratObj) <- clusteringResults
  figure <- DimPlot(seuratObj, label = label, cols = cols)
  return(figure)
}


