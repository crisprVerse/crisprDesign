#' Add a count-weighted specificity proxy score
#'
#' Computes a specificity proxy score for each gRNA based on the number of
#' genomic alignments with 0, 1, 2, or 3 mismatches. The score is intended as
#' a simple measure of off-target burden when empirical mismatch-tolerance
#' models are unavailable for the nuclease of interest.
#'
#' For each spacer, the score is calculated as:
#'
#' \deqn{
#' \mathrm{score} =
#' \frac{1}
#' {n_0 + w_1 n_1 + w_2 n_2 + w_3 n_3}
#' }
#'
#' where \eqn{n_0}, \eqn{n_1}, \eqn{n_2}, and \eqn{n_3} denote the number of
#' genomic alignments with 0, 1, 2, and 3 mismatches, respectively. If no
#' perfect match is identified (\eqn{n_0 = 0}), an additional penalty of 1 is
#' added to the denominator.
#'
#' Larger values indicate higher predicted specificity (lower off-target
#' burden), whereas smaller values indicate the presence of more closely
#' matching genomic sites. The score should be interpreted as a relative
#' ranking metric and not as a calibrated probability of off-target cleavage.
#'
#' @param guideSet A GuideSet object containing alignment information produced
#'   by \code{addSpacerAlignments()}.
#'
#' @param w1 Numeric weight assigned to alignments with one mismatch.
#'
#' @param w2 Numeric weight assigned to alignments with two mismatches.
#'
#' @param w3 Numeric weight assigned to alignments with three mismatches.
#'
#' @return
#' The input GuideSet with an additional column,
#' \code{score_specificity_proxy}, added to
#' \code{S4Vectors::mcols(guideSet)}.
#'
#' @details
#' This method provides a simple count-based approximation of guide
#' specificity. Unlike CFD or MIT specificity scores, it does not account for
#' mismatch position, mismatch identity, chromatin accessibility, or empirical
#' cleavage measurements. Consequently, it should be viewed as a heuristic
#' specificity proxy suitable for relative guide ranking.
#'
#' @seealso
#' \code{\link{addSpacerAlignments}}
#'
#'
#'
#' @export
addSpecificityProxyScore <- function(guideSet,
    w1=0.5,
    w2=0.05,
    w3=0.005
){

    aln <- alignments(guideSet)
    aln <- as.data.frame(S4Vectors::mcols(aln),
                         stringsAsFactors=FALSE)
    aln <- split(aln, f=aln$spacer)

    .getAggregateScore <- function(){
        vapply(aln, function(x){
            nmm <- x[["n_mismatches"]]
            nmm <- factor(nmm, levels=c(0,1,2,3))
            counts <- table(nmm)
            ws <- c(1,w1,w2,w3)
            total <- sum(counts*ws)

            if (counts["0"]==0){
                total <- total + 1
            } 
            return(1/total)
        }, FUN.VALUE=numeric(1))
    }
    
    scores <- .getAggregateScore()
    scores <- scores[match(as.character(spacers(guideSet)), names(scores))]
    S4Vectors::mcols(guideSet)[["score_specificity_proxy"]] <- scores
    return(guideSet)
}







