% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opsDesign.R
\name{getBarcodeDistanceMatrix}
\alias{getBarcodeDistanceMatrix}
\title{Get distance between query and target sets of barcodes}
\usage{
getBarcodeDistanceMatrix(
  queryBarcodes,
  targetBarcodes = NULL,
  binnarize = TRUE,
  min_dist_edit = NULL,
  dist_method = c("hamming", "levenshtein"),
  ignore_diagonal = TRUE,
  splitByChunks = FALSE,
  n_chunks = NULL
)
}
\arguments{
\item{queryBarcodes}{Character vector of DNA sequences or DNAStringSet.}

\item{targetBarcodes}{Optional character vector of DNA sequences 
or DNAStringSet. If NULL, distances will be calculated between
barcodes provided in \code{queryBarcodes}.}

\item{binnarize}{Should the distance matrix be made binnary?
TRUE by default. See details section.}

\item{min_dist_edit}{Integer specifying the minimum distance edit
required for barcodes to be considered dissimilar when
\code{binnarize=TRUE}, ignored otherwise.}

\item{dist_method}{String specifying distance method. Must be
either "hamming" (default) or "levenshtein".}

\item{ignore_diagonal}{When \code{targetBarcodes=NULL}, should the
diagonal distances be set to 0 to ignore self distances?
TRUE by default.}

\item{splitByChunks}{Should distances be calculated in a chunk-wise
manner? FALSE by default. Highly recommended when the set of query
barcodes is large to reduce memory footprint.}

\item{n_chunks}{Integer specifying the number of chunks to be used
when \code{splitByChunks=TRUE}. If NULL (default), number of chunks
will be chosen automatically.}
}
\value{
A sparse matrix of class \code{dgCMatrix} or \code{dsCMatrix}
    in which rows correspond to \code{queryBarcodes} and columns
    correspond to \code{targetBarcodes}. If \code{binnarize=TRUE},
    a value of 0 indicates that two barcodes have a distance 
    greater of equal to \code{min_dist_edit}, otherwise the value 
    is 1. If If \code{binnarize=FALSE}, values represent
    the actual calculated distances between barcodes.
}
\description{
Get distance between query and target sets of barcodes
}
\examples{
data(guideSetExample, package="crisprDesign")
guideSetExample <- addOpsBarcodes(guideSetExample)
barcodes <- as.character(guideSetExample$opsBarcode)
dist <- getBarcodeDistanceMatrix(barcodes, min_dist_edit=2)

}
\author{
Jean-Philippe Fortin
}
