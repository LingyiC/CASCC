% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmeansCenterModify.R
\name{kmeansCenterModify}
\alias{kmeansCenterModify}
\title{Kmeans with centers}
\usage{
kmeansCenterModify(res)
}
\arguments{
\item{res}{The CASCC object}
}
\description{
Kmeans with centers
}
\examples{
library("CASCC")
library("cafr")
library("Seurat")
data("Data_PDAC_peng_2k") 
data <- Data_PDAC_peng_2k
# Step 0 - Step 3
fs.res <- CASCC.featureSelection(data, inputDataType = "well-normalized", attr.raw = NULL, exponent.max = 10, exponent.min = 2, 
                                 generalSeedList = NULL, mc.cores = 1, topDEGs = 10, topAttr = 50,
                                 removeMT.RP.ERCC = TRUE, removeNonProtein = FALSE,
                                 topN.DEG.as.seed = 1,
                                 overlapN = 10)
adata  <- fs.res$adata; 
SeuratVersionCheck = as.numeric(strsplit(as.character(packageVersion("Seurat")), "\\\\.")[[1]][1])
if (SeuratVersionCheck == 4) {
  data <- adata[["RNA"]]@data
}else if (SeuratVersionCheck == 5) {
  data <- adata[["RNA"]]$data
}
data <- as.matrix(data)

# Step 4
res.predictK <- predictClusterK(data, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = FALSE)

# results list
res <- list()
res$fs.res <- fs.res
res$res.predictK <- res.predictK

# Step 5 
res$mainType.output <- kmeansCenterModify(res)

}
