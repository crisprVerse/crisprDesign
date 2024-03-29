% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/addOffTargetScores.R
\name{addOffTargetScores}
\alias{addOffTargetScores}
\alias{addOffTargetScores,GuideSet-method}
\alias{addOffTargetScores,PairedGuideSet-method}
\alias{addOffTargetScores,NULL-method}
\title{Add CFD and MIT scores to a \linkS4class{GuideSet} object.}
\usage{
addOffTargetScores(object, ...)

\S4method{addOffTargetScores}{GuideSet}(object, max_mm = 2, includeDistance = TRUE, offset = 0)

\S4method{addOffTargetScores}{PairedGuideSet}(object, max_mm = 2, includeDistance = TRUE, offset = 0)

\S4method{addOffTargetScores}{NULL}(object)
}
\arguments{
\item{object}{A \linkS4class{GuideSet} object or a 
\linkS4class{PairedGuideSet} object.
\code{crisprNuclease(object)} must be either using SpCas9 or CasRx.}

\item{...}{Additional arguments, currently ignored.}

\item{max_mm}{The maximimum number of mismatches between the spacer sequence
and the protospacer off-target sequence to be considered in the 
off-target score calculations. Off-targets with a number of 
mismatches greater than \code{max_mm} will be excluded; this is useful
if one wants to avoid the aggregated off-target scores to be driven by a
large number of off-targets that have low probability of cutting.}

\item{includeDistance}{Should a distance penalty for the MIT score be
included? TRUE by default.}

\item{offset}{Numeric value specifying an offset to add to the denominator
when calcuting the aggregated score (inverse summation formula).
0 by default.}
}
\value{
A \code{GuideSet} or a \code{PairedGuideSet} object with added 
    scores. The alignments annotation returned by \code{alignments(object)}
    will have additional column storing off-target scores. Those scores
    representing the off-target score for each gRNA and off-target pair.
    For SpCas9, a column containing an aggregated specificity off-target 
    score for each scoring method is added to the metadata columns 
    obtained by \code{mcols(object)}.
}
\description{
Add CFD and MIT off-target scores to a 
    \linkS4class{GuideSet} object. 
    Both the CFD and MIT methods are available for the SpCas9 nuclease.
    The CFD method is also available for the CasRx nuclease.
    Other nucleases are currently not supported.
}
\details{
See the \pkg{crisprScore} package for a description of the 
    different off-target scoring methods.
}
\examples{

data(guideSetExampleWithAlignments, package="crisprDesign")
gs <- guideSetExampleWithAlignments
gs <- addOffTargetScores(gs)

}
\seealso{
\code{link{addOnTargetScores}} to add on-target scores.
}
\author{
Jean-Philippe Fortin, Luke Hoberecht
}
