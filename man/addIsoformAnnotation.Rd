% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/addDistanceToTss.R,
%   R/addIsoformAnnotation.R
\name{addIsoformAnnotation}
\alias{addIsoformAnnotation}
\alias{addDistanceToTss,NULL-method}
\alias{addIsoformAnnotation,GuideSet-method}
\alias{addIsoformAnnotation,PairedGuideSet-method}
\alias{addIsoformAnnotation,NULL-method}
\title{Add isoform-specific annotation to a \linkS4class{GuideSet} object}
\usage{
addIsoformAnnotation(object, ...)

\S4method{addDistanceToTss}{NULL}(object)

\S4method{addIsoformAnnotation}{GuideSet}(object, tx_id)

\S4method{addIsoformAnnotation}{PairedGuideSet}(object, tx_id)

\S4method{addIsoformAnnotation}{NULL}(object)
}
\arguments{
\item{object}{A \linkS4class{GuideSet} object or a 
\linkS4class{PairedGuideSet} object.}

\item{...}{Additional arguments, currently ignored.}

\item{tx_id}{String specifiying Ensembl ID for the 
isoform transcript of interested. E.g. "ENST00000311936".}
}
\value{
A A \linkS4class{GuideSet} object or a 
    \linkS4class{PairedGuideSet} object.with the 
    following added columns: \code{percentCDS},
    \code{percentCodingIsoforms}, and 
    \code{isCommonCodingExon}. The column values are
    specific to the transcript specified by \code{tx_id}.
    The \code{percentCDS} column indicates at what percentage
    of the coding sequence the gRNA is cutting. The
    column \code{percentCodingIsoforms} indicates the 
    percentage of coding isoforms that are targeted
    by the gRNA. The column \code{isCommonCodingExon}
    indicates whether or not the exon targetd by the
    gRNA is common to all isoforms for the gene.
}
\description{
Add isoform-specific annotation to a
     \linkS4class{GuideSet} object.
}
\examples{
data(guideSetExampleFullAnnotation)
tx_id <- "ENST00000538872"
gs <- guideSetExampleFullAnnotation
gs <- addIsoformAnnotation(gs, tx_id)

}
\author{
Jean-Philippe Fortin
}
