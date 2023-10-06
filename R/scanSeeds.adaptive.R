##' Scan seeds in adaptive manner
##' 
##' 
##' 
##' @param data A matrix of gene expression data. Rows are genes and columns
##' are cells.
##' @param markers A list of genes that are used as seeds.
##' @param epsilon Threshold of convergence. By default, epsilon = 1e-4.
##' @param verbose When TRUE, it will show the top 20 genes of the metagene in
##' each iteration.
##' @param exponent.max max exponent value.
##' @param exponent.min min exponent value.
##' @param seedTopN The threshold of divergence.
##' @param compareNth Maximize the strength of the compareNth-ranked genes in
##' the converged attractor.
##' @param mc.cores Number of cores to be used. Default is 1.
##' @return
##' @author Lingyi Cai
##' @seealso \code{\link{findAttractor.adaptive}}
##' @references
##' @examples
##' 
##' 
##' @export
 
scanSeeds.adaptive <- function(data, markers, epsilon = 1e-4, verbose = FALSE, exponent.max = 10, exponent.min = 2, seedTopN = min(round(nrow(data)*0.005), 50), compareNth = 10, mc.cores = NULL){

  if(is.null(mc.cores) == TRUE){
    library("cafr")
    data <- as.matrix(data)
    seedList <- markers
    attractorObjects <- list()
    attractorTouchedGenes <- list()
    history_exponent <- list()
    final.exponent <- list()
    i = 0
    # genesOrder <- rownames(data)

    scannedSeeds <- c()
    attractorSignatureList <- list()
    runSeed <- list()
    while(identical(seedList, character(0)) == FALSE){
      cat(length(seedList),"seeds to be scanned.\r")
      flush.console()
      seed = seedList[1]
      i = i + 1
      runSeed[[i]] <- seed
      tmp <- NULL
      tryCatch({
        tmp <- findAttractor.adaptive(data, seed, exponent.max = exponent.max, exponent.min = exponent.min, step.large = 1, step.small = 0.1, dominantThreshold = 0.2, seedTopN = seedTopN, compareNth = compareNth, verbose = FALSE, epsilon = epsilon, minimize = FALSE, first.idx = 1, second.idx = 2)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

      if (is.null(tmp)){
        tmp <- list()
        tmp$attr <- rep(0, nrow(data))
        names(tmp$attr ) <- rownames(data)
      }
      attractorObjects[[i]] <- tmp$attr
      attractorTouchedGenes[[i]] <- tmp$scannedSeeds
      history_exponent[[i]] <- tmp$history_exponent
      final.exponent[[i]] <- tmp$final.exponent
      # seedList <- setdiff(seedList, c(names(attractorObjects[[i]]), seed))
      # scannedSeeds <- c(scannedSeeds, c(names(attractorObjects[[i]]), seed))
      seedList <- setdiff(seedList, tmp$scannedSeeds)
      scannedSeeds <- c(scannedSeeds, tmp$scannedSeeds)
      attractorSignatureList[[i]] <- attractorObjects[[i]]
    }
    res <- list(attractorSignatureList=attractorSignatureList, attractorTouchedGenes = attractorTouchedGenes, runSeed = runSeed, history_exponent=history_exponent, final.exponent=final.exponent)
  }else{
    library("cafr")
    data <- as.matrix(data)
    seedList <- markers
    attractorObjects <- list()
    attractorTouchedGenes <- list()
    history_exponent <- list()
    final.exponent <- list()


    scannedSeeds <- c()
    attractorSignatureList <- list()
    runSeed <- list()
    ## parallel
    parallel.res <- parallel::mclapply(
      seedList,
      function(seed) {
        tryCatch(
          findAttractor.adaptive(data, seed, exponent.max = exponent.max, exponent.min = exponent.min, step.large = 1, step.small = 0.1,
                                 dominantThreshold = 0.2, seedTopN = seedTopN, compareNth = compareNth, verbose = FALSE, epsilon = epsilon,
                                 minimize = FALSE, first.idx = 1, second.idx = 2),
          error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      },
      mc.cores = mc.cores
    )
    for (i in seq(1, length(parallel.res))) {
      tmp <- parallel.res[[i]]
      if (is.null(tmp)){
        tmp <- list()
        tmp$attr <- rep(0, nrow(data))
        names(tmp$attr ) <- rownames(data)
      }
      attractorObjects[[i]] <- tmp$attr
      attractorTouchedGenes[[i]] <- tmp$scannedSeeds
      history_exponent[[i]] <- tmp$history_exponent
      final.exponent[[i]] <- tmp$final.exponent
      # seedList <- setdiff(seedList, c(names(attractorObjects[[i]]), seed))
      # scannedSeeds <- c(scannedSeeds, c(names(attractorObjects[[i]]), seed))
      seedList <- setdiff(seedList, tmp$scannedSeeds)
      scannedSeeds <- c(scannedSeeds, tmp$scannedSeeds)
      attractorSignatureList[[i]] <- attractorObjects[[i]]
    }
    res <- list(attractorSignatureList=attractorSignatureList, attractorTouchedGenes = attractorTouchedGenes, runSeed = c("parallel: ", markers), history_exponent=history_exponent, final.exponent=final.exponent)


    ## make the parallel resulst be consistent with single core
    attrs <- res$attractorSignatureList
    target <- res$attractorTouchedGenes
    seeds <- res$runSeed[-1]
    removed_index <- c()
    for (i in 2:length(target)) {
      if(seeds[i] %in% unlist(target[1:(i-1)])){
        removed_index <- c(removed_index, i)
      }
    }

    if(length(removed_index) > 0){
      target[removed_index] <- NULL
      attrs[removed_index] <- NULL
      res$attractorSignatureList_parallel_old <- res$attractorSignatureList
      res$attractorSignatureList_parallel_removed_index <- removed_index
      res$attractorSignatureList <- attrs
    }
  }



  return(res)
}




findAttractor.adaptive_forMultiCores <- function(data, seed, epsilon = 1e-4, verbose = FALSE, exponent.max = 5, exponent.min = 2){
      tmp <- NULL
      tryCatch({
        tmp <- findAttractor.adaptive(data, seed, exponent.max = exponent.max, exponent.min = exponent.min, step.large = 1, step.small = 0.1, dominantThreshold = 0.2, seedTopN = 100, compareNth = compareNth, verbose = FALSE, epsilon = epsilon, minimize = FALSE, first.idx = 1, second.idx = 2)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

      if (is.null(tmp)){
        tmp <- list()
        tmp$attr <- rep(0, nrow(data))
        names(tmp$attr ) <- rownames(data)
      }

      return(tmp$attr)
}

# load("/drive_T7/2021Research/Peng/data/adata_byPatient/adata.N3.RData")
# data <- adata[["RNA"]]@data
# attractorSignatureList <- mclapply(seedList,function(seed, data){findAttractor.adaptive_forMultiCores(data,seed)}, data = data, mc.cores = 2)
# res <- list(attractorSignatureList=attractorSignatureList)

