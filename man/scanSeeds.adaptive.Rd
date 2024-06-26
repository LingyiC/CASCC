% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanSeeds.adaptive.R
\name{scanSeeds.adaptive}
\alias{scanSeeds.adaptive}
\title{Scan seeds in adaptive manner}
\usage{
scanSeeds.adaptive(
  data,
  markers,
  epsilon = 1e-04,
  verbose = FALSE,
  exponent.max = 10,
  exponent.min = 2,
  seedTopN = min(round(nrow(data) * 0.005), 50),
  compareNth = 10,
  mc.cores = NULL
)
}
\arguments{
\item{data}{A matrix of gene expression data. Rows are genes and columns
are cells.}

\item{markers}{A list of genes that are used as seeds.}

\item{epsilon}{Threshold of convergence. By default, epsilon = 1e-4.}

\item{verbose}{When TRUE, it will show the top 20 genes of the metagene in
each iteration.}

\item{exponent.max}{max exponent value.}

\item{exponent.min}{min exponent value.}

\item{seedTopN}{The threshold of divergence.}

\item{compareNth}{Maximize the strength of the compareNth-ranked genes in
the converged attractor.}

\item{mc.cores}{Number of cores to be used. Default is 1.}
}
\value{
\item{attractorSignatureList}{Attractors.} \item{attractorTouchedGenes}{Skipped 
seeds.} \item{runSeed}{Scanned seeds.} 
\item{history_exponent}{History
exponent.} \item{final.exponent}{Final exponent.}
}
\description{
Scan seeds in adaptive manner
}
\examples{
library("CASCC")
library("cafr")
data("Data_PDAC_peng_2k") 
# Find LUM-seeded and RGS5-seeded attractor. LUM gene is a marker of fibroblast cells. RGS5 gene is a marker of pericyte cells.
attr <- scanSeeds.adaptive(Data_PDAC_peng_2k, c("LUM", "RGS5"))
}
\seealso{
\code{\link{findAttractor.adaptive}}
}
\author{
Lingyi Cai
}
