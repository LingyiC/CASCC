# CASCC

CASCC, a clusteing method desined to improve biological accuracy using gene co-expression features identified using an unsupervised adaptive attractor algorithm. 

The users can refer to our paper xxxx for a detailed description of the modeling and applications.

- [Overview](#Overview)
- [Tutorials](#Tutorials)
    - [Installation](#Installation)
    - [Quick start](#Quick-start)
    - [Run CASCC step by step](#Run-CASCC-step-by-step)
- [Other files](#Other-files)
    - [Toy data](#Toy-data)
- [Session info](#Session-info) 

# Tutorials
## Installation

You can install CASCC from GitHub with:

```R
require(devtools)
install_github("LingyiC/CASCC")
install_github("weiyi-bitw/cafr")  # The 'cafr' R package is required for running CASCC.
```
## Quick-start
```R
library(CASCC)
## load the toy data
data("Data_Tissues_li_2k") 

## run CASCC
res <- CASCC::run.CASCC(data) # inputDataType `data` is normalized and log-transformed. 
# groundTruth = NULL, attr.raw = NULL, inputDataType = "well-normalized"
# If users already have attractor scanning results, they can input `attr.raw` to save time.

## plot the results
adata <- res$res.predictK$adata # feature-selected UMAP
labels <- res$mainType.output$clusteringResults # clustering output: k-means with attractor centers
Idents(adata) <- labels
DimPlot(adata)

## results list
attr.raw <- res$fs.res$attr.raw # raw attractor scanning results, can be used if rerun CASCC
attrs <- res$fs.res$attrs # attractors removed identical ones
finalAttrs <- res$fs.res$finalAttrs # filtered attractors 

```
## Run-CASCC-step-by-step
```R
library("CASCC")
# Step 0 - Step 3
fs.res <- CASCC.featureSelection(data, inputDataType = "well-normalized", attr.raw = NULL, exponent.max = 10, exponent.min = 2, 
                                 generalSeedList = NULL, mc.cores = 1, topDEGs = 10, topAttr = 50,
                                 removeMT.RP.ERCC = TRUE, removeNonProtein = FALSE,
                                 topN.DEG.as.seed = 1,
                                 overlapN = 10)
adata  <- fs.res$adata; data <- adata[["RNA"]]@data; data <- as.matrix(data)

# Step 4
res.predictK <- predictClusterK(data, fs.res, method = "ward.D2", index = "silhouette", min.nc_fix = FALSE)

# Step 5 - fixed centers, mainType.output clustering results: res$maintype.output$clusteringResults
tmp <- kmeansCenterModify(res)


# results list
res <- list()
res$fs.res <- fs.res
res$res.predictK <- res.predictK
res$mainType.output <- tmp
```


# Other-files

## Toy-data
Generate the toy data
```R
load("./Data_Tissues_li_full.RData") # download from https://hemberg-lab.github.io/scRNA.seq.datasets/,  data <- rds@assays[["data"]]@listData[["logcounts"]]
library("Seurat")
adata <- CreateSeuratObject(counts = data, verbose = F)
# adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000, verbose = F)

data <- adata[["RNA"]]@data[VariableFeatures(adata), ]
data <- as.matrix(data)

save(data, file = "./CASCC/data/Data_Tissues_li_2k.RData")

library(CASCC)
## load the toy data
data("Data_Tissues_li_2k")

## run CASCC
res <- CASCC::run.CASCC(data, groundTruth = NULL, min.nc_fix = FALSE, attr.raw = NULL)

## plot the results
adata <- res$res.predictK$adata # feature-selected UMAP
labels <- res$mainType.output$clusteringResults # clustering output: k-means with attractor centers
Idents(adata) <- labels
DimPlot(adata)

save(res, file = "./CASCC/data/Res_Tissues_li_2k.RData")
```

# Session-info
Session info to replicate the clustering results of Table S1 datasets
```
R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cafr_0.312    dplyr_1.1.2   aricode_1.0.2 stringr_1.5.0

loaded via a namespace (and not attached):
  [1] Seurat_4.3.0                Rtsne_0.16                 
  [3] colorspace_2.1-0            deldir_1.0-6               
  [5] ellipsis_0.3.2              ggridges_0.5.4             
  [7] XVector_0.30.0              GenomicRanges_1.42.0       
  [9] spatstat.data_3.0-0         leiden_0.4.3               
 [11] listenv_0.9.0               ggrepel_0.9.3              
 [13] fansi_1.0.4                 codetools_0.2-19           
 [15] splines_4.0.5               polyclip_1.10-4            
 [17] jsonlite_1.8.4              ica_1.0-3                  
 [19] cluster_2.1.4               png_0.1-8                  
 [21] uwot_0.1.14                 shiny_1.7.4                
 [23] sctransform_0.3.5           spatstat.sparse_3.0-0      
 [25] compiler_4.0.5              httr_1.4.5                 
 [27] SeuratObject_4.1.3          Matrix_1.5-3               
 [29] fastmap_1.1.1               lazyeval_0.2.2             
 [31] limma_3.46.0                cli_3.6.1                  
 [33] later_1.3.0                 htmltools_0.5.5            
 [35] tools_4.0.5                 igraph_1.4.2               
 [37] gtable_0.3.3                glue_1.6.2                 
 [39] GenomeInfoDbData_1.2.4      RANN_2.6.1                 
 [41] reshape2_1.4.4              Rcpp_1.0.10                
 [43] scattermore_0.8             Biobase_2.50.0             
 [45] vctrs_0.6.2                 nlme_3.1-162               
 [47] spatstat.explore_3.0-6      progressr_0.13.0           
 [49] lmtest_0.9-40               spatstat.random_3.1-3      
 [51] globals_0.16.2              mime_0.12                  
 [53] miniUI_0.1.1.1              lifecycle_1.0.3            
 [55] irlba_2.3.5.1               goftest_1.2-3              
 [57] future_1.32.0               zlibbioc_1.36.0            
 [59] MASS_7.3-58.3               zoo_1.8-12                 
 [61] scales_1.2.1                promises_1.2.0.1           
 [63] MatrixGenerics_1.2.1        spatstat.utils_3.0-1       
 [65] NbClust_3.0.1               parallel_4.0.5             
 [67] SummarizedExperiment_1.20.0 RColorBrewer_1.1-3         
 [69] reticulate_1.28             pbapply_1.7-0              
 [71] gridExtra_2.3               ggplot2_3.4.2              
 [73] stringi_1.7.12              S4Vectors_0.28.1           
 [75] BiocGenerics_0.36.1         GenomeInfoDb_1.26.7        
 [77] rlang_1.1.0                 pkgconfig_2.0.3            
 [79] matrixStats_0.63.0          bitops_1.0-7               
 [81] lattice_0.20-45             tensor_1.5                 
 [83] ROCR_1.0-11                 purrr_1.0.1                
 [85] patchwork_1.1.2             htmlwidgets_1.6.2          
 [87] cowplot_1.1.1               tidyselect_1.2.0           
 [89] parallelly_1.35.0           RcppAnnoy_0.0.20           
 [91] plyr_1.8.8                  magrittr_2.0.3             
 [93] R6_2.5.1                    IRanges_2.24.1             
 [95] generics_0.1.3              DelayedArray_0.16.3        
 [97] pillar_1.9.0                fitdistrplus_1.1-11        
 [99] abind_1.4-5                 survival_3.5-3             
[101] RCurl_1.98-1.10             sp_1.6-0                   
[103] tibble_3.2.1                future.apply_1.10.0        
[105] CASCC_0.1.0                 KernSmooth_2.23-20         
[107] utf8_1.2.3                  spatstat.geom_3.0-6        
[109] plotly_4.10.1               grid_4.0.5                 
[111] data.table_1.14.8           digest_0.6.31              
[113] xtable_1.8-4                tidyr_1.3.0                
[115] httpuv_1.6.9                stats4_4.0.5               
[117] munsell_0.5.0               viridisLite_0.4.1          
```