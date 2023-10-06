## orginal file /drive_T7/2021Research/CASCC_0317/CASCC_0317/CASCC_LC/attractorLib2.R


##' Collapse AttractorList
##' 
##' This function helps remove identical or similar attractors from a list of
##' attractors.
##' 
##' 
##' @param attractorList A list of attractors.
##' @param NumTopFeature Number of top features to be used in the attractors'
##' comparison.
##' @param verbose Whether to print out the progress.
##' @param overlapN The number of top features that are allowed to overlap
##' between two attractors.
##' @return
##' @author Lingyi Cai
##' @seealso
##' @references
##' @examples
##' 
##' 
##' 
##' @export



collapseAttractorList <- function(attractorList, NumTopFeature=50, outputLength=NULL, verbose=TRUE, overlapN = 50, removeDominant = FALSE){
## LC - remove the similar attractors
	attractorCount <-0
	allSignatureList<-list()

	if(verbose == TRUE)
		cat("Collapsing the found attractors...")

	if(length(attractorList) > 0)
	if(length(attractorList$attractorSignatureList) > 0){
		for(j in 1:length(attractorList$attractorSignatureList)){

			if(is.null(outputLength) == TRUE){
				pendingAttractor<-sort( attractorList$attractorSignatureList[[j]], decreasing=TRUE )
			}


			doesAttractorExist<-FALSE
			if(length(allSignatureList) > 0){
				for(k in 1:length(allSignatureList)){
					doesAttractorExist<-areVectorsSimilar(names(pendingAttractor)[1:NumTopFeature], names(allSignatureList[[k]])[1:NumTopFeature], overlapN = overlapN)
					if(doesAttractorExist == TRUE){
						## 110922, LC added the comparison
						isMoreStable <- (pendingAttractor[1] - pendingAttractor[2]) < (allSignatureList[[k]][1] - allSignatureList[[k]][2] )
						if(isMoreStable == TRUE){
						allSignatureList[[k]] <- pendingAttractor
						}
						break
					}

				}
			}

			if(doesAttractorExist == FALSE){
				attractorCount=attractorCount+1
				allSignatureList[[attractorCount]]<-pendingAttractor
			}
		}
	}

	if(verbose == TRUE)
		cat("DONE.\n")

	allSignatureList[which((unlist(lapply(allSignatureList, sum))) %in% c(0, NumTopFeature))] <- NULL ## remove 111 or 000 attractor
	if(removeDominant == TRUE){
		allSignatureList[which((unlist(lapply(allSignatureList, function(x) {round(as.numeric(x[1]), 3)}))) == 1 )] <- NULL ## remove attractor with first gene = 1, dominant
	}

	return(allSignatureList)
}


attractorModify <- function(attrs, adata, expressionLevel = 0, percentage = 0.02, overlapN = 10){

  attr.last.num = 0
  while (attr.last.num != length(attrs)) {
    test<- list(); test$attractorSignatureList <- attrs
    attrs <- collapseAttractorList(test, NumTopFeature=50, overlapN = overlapN, verbose = F) # original = 10
    attr.last.num <- length(attrs)
  }

  ## remove attractors do not express in at least 2% cells, by default we skip this step 
  top1_genes <- unlist(lapply(attrs, function(x){names(x[1])}))
  tmp <- which(sapply(top1_genes, function(x){length(which(adata[["RNA"]]@data[x, ] > expressionLevel))}) < ncol(adata)*percentage)
  attrs[tmp] <- NULL
  # tmp <- which(lapply(attrs, function(x){x[1]}) > 0.99)
  attrs[tmp] <- NULL
  return(attrs)
}






#Collapse
areVectorsSimilar <-function(vectorA, vectorB, verbose=FALSE, overlapN=5){
	isSimilar=FALSE

	if(length(intersect(vectorA,vectorB))>=overlapN)
		isSimilar=TRUE

	if(verbose==TRUE){
		cat(vectorA[1],vectorB[1], length(intersect(vectorA,vectorB)) ,"\n")
	}
	return(isSimilar)
}
