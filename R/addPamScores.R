#' @title Add PAM scores to a \linkS4class{GuideSet} object.
#' @description Add PAM scores to a \linkS4class{GuideSet} object
#'     based on the \linkS4class{CrisprNuclease} object stored in
#'     the \linkS4class{GuideSet} object. PAM scores indicate nuclease
#'     affinity (recognition) to different PAM sequences.
#'     A score of 1 indicates a PAM sequence that is fully
#'     recognized by the nuclease.
#' 
#' @param guideSet A \linkS4class{GuideSet} or a 
#'     \linkS4class{PairedGuideSet} object.
#' 
#' @return \code{guideSet} with an appended \code{score_pam} column in
#'     \code{mcols(guideSet)}.
#'
#' @examples
#' 
#' # Using character vector as input:
#' data(enAsCas12a, package="crisprBase")
#' gs <- findSpacers("CCAACATAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTAGCCGGGTGTTA",
#'                   canonical=FALSE,
#'                   crisprNuclease=enAsCas12a)
#' gs <- addPamScores(gs)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom crisprBase weights
#' @importFrom S4Vectors mcols<-
addPamScores <- function(guideSet){
    guideSet <- .validateGuideSetOrPairedGuideSet(guideSet)
    if (.isGuideSet(guideSet)){
        out <- addPamScores_guideset(guideSet)
    } else if (.isPairedGuideSet(guideSet)){
        unifiedGuideSet <- .pairedGuideSet2GuideSet(guideSet)
        unifiedGuideSet <- addPamScores_guideset(unifiedGuideSet)
        out <- .addColumnsFromUnifiedGuideSet(guideSet,
                                              unifiedGuideSet)
    }
    return(out)
}


#' @importFrom crisprBase weights
#' @importFrom S4Vectors mcols<-
addPamScores_guideset <- function(guideSet){
    guideSet <- .validateGuideSet(guideSet)
    crisprNuclease <- crisprNuclease(guideSet)
    pamScoreWeights <- crisprBase::weights(crisprNuclease, expand=TRUE)
    matchingIndices <- match(pams(guideSet, as.character=TRUE),
                             names(pamScoreWeights))
    scores <- round(pamScoreWeights[matchingIndices], 2)
    S4Vectors::mcols(guideSet)[["score_pam"]] <- scores
    return(guideSet)
}
