##' Run CASCC in one line of code
##'
##' Run the CASCC package with one line of code.
##'
##'
##' @param data.lognormcount Input gene expression matrix.
##' @param groundTruth.K Default \code{NULL}. The number of clusters. By
##' default, CASCC will estimate the number of clusters. Users can manually set
##' the number of clusters if they have prior knowledge about the number of clusters.
##' @param attr.raw Default \code{NULL}. If users already run CASCC and have attractor
##' scanning results \code{CASCC.output$fs.res$attrs.raw}, they can provide the
##' list of attractors.
##' @param inputDataType Input data type. Default is "well-normalized".
##' Options: "well-normalized", "needLog2" or "UMI". If the input matrix is
##' normalized and log-transformed, use "well-normalized". If the input matrix
##' is normalized but not log-transformed, use "needLog2". If the input matrix
##' is UMI-count, use "UMI".
##' @param mc.cores Default \code{NULL}. The number of cores to be used. By
##' default, CASCC use 2 cores.
##' @param removeMT.RP.ERCC Default \code{TRUE}. By default, mitochondrial genes,
##' ribosomal genes, and spike-in genes will be removed.
##' @param removeNonProtein Default \code{FALSE}. By default, non-protein coding
##' genes will not be removed.To facilitate the analysis, users can remove them.
##' @param topN.DEG.as.seed Default \code{1}. By default, the top 1 DEG of each
##' initial cluster will be used as seeds in the following analysis.
##' @param generalSeedList Seed genes set by users. Default is NULL. For example,
##' \code{generalSeedList = c("AMBP", "KRT19", "PRSS1", "CHGB", "RGS5", "LUM", "CDH5", "AIF1", "CD3D", "MS4A1")}.
##' @param exponent.max Default \code{10}.To facilitate the analysis of large dataset,
##' users can set a smaller value, such as 5.
##' @param exponent.min Default \code{2}.
##' @param min.diff.pct Default \code{0.5}. Seurat - Only test genes that show a minimum difference 
##' in the fraction of detection between the two groups. Set to 0.5 by default.
##' @param min.pct Default \code{0.25}.  Seurat - Only test genes that are detected 
##' in a minimum fraction of min.pct cells in either of the two populations.
##' Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.25.
##' @param logfc.threshold Default \code{0.25}. Seurat - The log fold change threshold. 
##' @param topDEGs Default \code{10}. By default, the top 10 DEGs of each
##' initail cluster will be used as features in the following analysis.
##' @param topAttr Default \code{50}. By default, the top 50 genes of each
##' attractor will be used as features in the following analysis.
##' @param min.nc_fix Default \code{FALSE}. If \code{TRUE}, the minimum number
##' of clusters will be fixed to be 3. If \code{FALSE}, the minimum number of
##' clusters will be suggested by CASCC attractor output.
##' @param if_PCA_input Default \code{TRUE}. If \code{TRUE}, the input will be
##' PCA transformed. If \code{FALSE}, the input will be reduced to 2D UMAP
##' space. If \code{"feature-matrix"}, the input will be the feature matrix.
##' @param nPCA Default \code{30}. The number of PCA components to be used.
##' @param overlapN Default \code{10}. The number of overlapping genes between
##' different attractors to define duplicated attractors.
##' @param attractorModify.percentage Default \code{0}.
##' Remove attractor expressed in less than \code{attractorModify.percentage} percent of cells. This is disabled by default.
##' @return
##' \item{fs.res}{Output of \code{\link{CASCC.featureSelection}}. Step
##' 1 to step 3 of CASCC.}
##' \item{res.predictK}{Output of
##' \code{\link{predictClusterK}}. Step 4 of CASCC.}
##' \item{mainType.output}{Output of \code{\link{kmeansCenterModify}}. Step 5
##' of CASCC. Kmeans with fixed centers. }
##' \item{res.kmeans}{Output
##' of run.kmeans function. Kmeans with random centers.
##' }
##'
##' @examples
##' # load the toy dataset
##' data("Data_PDAC_peng_2k")
##' dim(Data_PDAC_peng_2k) # 400 cells and 2000 genes
##'
##' # run CASCC
##' res <- CASCC::run.CASCC(Data_PDAC_peng_2k) # inputDataType `data` is normalized and log-transformed.
##'
##' # clustering results
##' labels <- res$mainType.output$clusteringResults
##'
##' @export

run.CASCC <- function(data.lognormcount, groundTruth.K = NULL, attr.raw = NULL, inputDataType = "well-normalized", mc.cores = 2,
                      removeMT.RP.ERCC = TRUE, removeNonProtein = FALSE,
                      min.diff.pct = 0.5, min.pct = 0.25, logfc.threshold = 0.25,
                      topN.DEG.as.seed = 1, generalSeedList = NULL, exponent.max = 10, exponent.min = 2, topDEGs = 10, topAttr = 50,
                      min.nc_fix = FALSE, if_PCA_input = TRUE, nPCA = 30,
                      overlapN = 10, attractorModify.percentage = 0){

    SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])
    set.seed(1)


    res <- list()

    # Step 0 - Step 3
    data.lognormcount <- as.matrix(data.lognormcount)
    fs.res <- CASCC.featureSelection(data.lognormcount, inputDataType = inputDataType, attr.raw = attr.raw, exponent.max = exponent.max, exponent.min = exponent.min, generalSeedList = generalSeedList, mc.cores = mc.cores, topDEGs = topDEGs, topAttr = topAttr,
                                     removeMT.RP.ERCC = removeMT.RP.ERCC, removeNonProtein = removeNonProtein,
                                     min.diff.pct = min.diff.pct, min.pct = min.pct, logfc.threshold = logfc.threshold,
                                     topN.DEG.as.seed = topN.DEG.as.seed,
                                     overlapN = overlapN, attractorModify.percentage = attractorModify.percentage)
    adata  <- fs.res$adata;


    if (SeuratVersionCheck == 4) {
      data.lognormcount <- adata[["RNA"]]@data
    }else if (SeuratVersionCheck == 5) {
      data.lognormcount <- adata[["RNA"]]$data
    }
    data.lognormcount <- as.matrix(data.lognormcount)

    # Step 4
    cat("Step 4: estimating the number of clusters...", "\n")
    res.predictK <- predictClusterK(data.lognormcount, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = min.nc_fix, if_PCA_input = if_PCA_input, nPCA = nPCA)

    # Step 5 - random centers
    if(is.null(groundTruth.K) == TRUE){
      res.kmeans <- run.kmeans(data.lognormcount[fs.res$features, ], res.predictK$clusterK)
    }else{
      res.kmeans <- run.kmeans(data.lognormcount[fs.res$features, ], groundTruth.K)
    }
    res$fs.res <- fs.res
    res$res.predictK <- res.predictK
    res$res.kmeans <- res.kmeans # kmeans without centers

    # Step 5 - fixed centers, mainType.output clustering results: res$maintype.output$clusteringResults
    cat("Step 5: performing clustering...", "\n")
    res$mainType.output <- kmeansCenterModify(res)

    cat("Done.", "\n")
    return(res)
}
