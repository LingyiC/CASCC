##' Run CASCC in one line of code
##'
##' Run the CASCC package with one line of code.
##'
##'
##' @param data.lognormcount Input gene expression matrix.
##' @param groundTruth.K Default \code{NULL}. The number of clusters. By
##' default, CASCC will estimate the number of clusters. Users can manually set
##' the number of clusters if they prefer.
##' @param attr.raw Default \code{NULL}. If users already have attractor
##' scanning results \code{CASCC.output$fs.res$attrs.raw}, they can provide the
##' list of attractors.
##' @param inputDataType Default \code{"well-normalized"}.
##' @param min.nc_fix Default \code{FALSE}. If \code{TRUE}, the minimum number
##' of clusters will be fixed to be 3. If \code{FALSE}, the minimum number of
##' clusters will be suggested by CASCC attractor output.
##' @param exponent.max Default \code{10}.
##' @param exponent.min Default \code{2}.
##' @param mc.cores Default \code{NULL}. The number of cores to be used. By
##' default, CASCC use only 1 core.
##' @param topDEGs Default \code{10}. By default, the top 10 DEGs of each
##' cluster will be used as features in the following analysis.
##' @param topAttr Default \code{50}. By default, the top 50 genes of each
##' attractor will be used as features in the following analysis.
##' @return \item{fs.res}{Output of \code{\link{CASCC.featureSelection}}. Step
##' 1 to step 3 of CASCC.} \item{res.predictK}{Output of
##' \code{\link{predictClusterK}}. Step 4 of CASCC.} \item{res.kmeans}{Output
##' of run.kmeans function. Kmeans with random centers.
##' }
##' \item{mainType.output}{Output of \code{\link{kmeansCenterModify}}. Step 5
##' of CASCC. Kmeans with fixed centers. }
##'
##' @references
##' @examples
##'
##'
##' @export

run.CASCC <- function(data.lognormcount, groundTruth.K = NULL, attr.raw = NULL, inputDataType = "well-normalized", mc.cores = 2,
                      removeMT.RP.ERCC = TRUE, removeNonProtein = FALSE,
                      topN.DEG.as.seed = 1, generalSeedList = NULL, exponent.max = 10, exponent.min = 2, topDEGs = 10, topAttr = 50,
                      min.nc_fix = FALSE,
                      overlapN = 10, attractorModify.percentage = 0){
    set.seed(1)
    res <- list()

    # Step 0 - Step 3
    fs.res <- CASCC.featureSelection(data.lognormcount, inputDataType = inputDataType, attr.raw = attr.raw, exponent.max = exponent.max, exponent.min = exponent.min, generalSeedList = generalSeedList, mc.cores = mc.cores, topDEGs = topDEGs, topAttr = topAttr,
                                     removeMT.RP.ERCC = removeMT.RP.ERCC, removeNonProtein = removeNonProtein,
                                     topN.DEG.as.seed = topN.DEG.as.seed,
                                     overlapN = overlapN, attractorModify.percentage = attractorModify.percentage)
    adata  <- fs.res$adata; data.lognormcount <- adata[["RNA"]]@data; data.lognormcount <- as.matrix(data.lognormcount)

    # Step 4
    res.predictK <- predictClusterK(data.lognormcount, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = min.nc_fix)

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
    res$mainType.output <- kmeansCenterModify(res)


    return(res)
}
