##' Filter genes
##'
##' This function helps remove genes which are not relevant for characterizing
##' cell types.
##'
##'
##' @param data Normalized data matrix.
##' @param min.cells default is 0. The genes that are not expressed in any
##' cells are removed.
##' @param removeNonProtein When \code{TRUE}, only keep protein coding genes.
##' @param removeMT.RP.ERCC When \code{TRUE}, filter out mitochondrial genes,
##' ribosomal genes, and spike-in gene.
##' @note First, the rows (features) with less than (min.cells + 1) cell that
##' have an expression level > 0 are excluded. Second, to remove genes not 
##' relevant for characterizing cell types, CASCC filtered out mitochondrial genes, 
##' ribosomal genes, and spike-in genes from the analysis by default.
##' @author Lingyi Cai
##' @seealso \code{\link{CASCC.featureSelection}}
##' @references
##' @examples
##' library("CASCC")
##' data("Data_PDAC_peng_2k") 
##' data <- filterGenes(Data_PDAC_peng_2k)
##' @export
filterGenes <- function(data, min.cells = 0.05*ncol(data), removeNonProtein = FALSE, removeMT.RP.ERCC = FALSE){
  # min.cells = 0.05*ncol(data)
  ## Include features detected in at least this many cells.
  filter.expr <- data[Matrix::rowSums(data > 0) > min.cells, ]

  ## remove protein
  if (removeNonProtein == TRUE){
    map <- data.frame(annotables::grch38[, c("ensgene", "symbol", "biotype")])
    nonProtein <- map[which(!map$biotype %in% "protein_coding"), ]$symbol
    Protein <- map[which(map$biotype %in% "protein_coding"), ]$symbol
    nonProtein <- unique(setdiff(nonProtein, Protein)); rm(Protein)
    genes <- setdiff(rownames(filter.expr), nonProtein)
    filter.expr <- filter.expr[genes, ]
  }

  if (removeMT.RP.ERCC == TRUE){
  ## remove mitochondrial genes, ribosomal genes, spike-in genes
  filter.expr <- filter.expr[!grepl(pattern = "^MT-|^RPS|^RPL|^ERCC-|^Rps|^Rpl", rownames(filter.expr)), ]
  }


  return(filter.expr)
}
