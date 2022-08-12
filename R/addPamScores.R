#' @title Add PAM scores to a \linkS4class{GuideSet} object.
#' @description Add PAM scores to a \linkS4class{GuideSet} object
#'     based on the \linkS4class{CrisprNuclease} object stored in
#'     the \linkS4class{GuideSet} object. PAM scores indicate nuclease
#'     affinity (recognition) to different PAM sequences.
#'     A score of 1 indicates a PAM sequence that is fully
#'     recognized by the nuclease.
#' 
#' @param object A \linkS4class{GuideSet} or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param ... Additional arguments, currently ignored.
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
#' @export
#' @importFrom crisprBase weights
#' @importFrom S4Vectors mcols<-
#' @rdname addPamScores
setMethod("addPamScores", "GuideSet", function(object){
    object <- .validateGuideSet(object)
    crisprNuclease <- crisprNuclease(object)
    pamScoreWeights <- crisprBase::weights(crisprNuclease, expand=TRUE)
    matchingIndices <- match(pams(object, as.character=TRUE),
                             names(pamScoreWeights))
    scores <- round(pamScoreWeights[matchingIndices], 2)
    S4Vectors::mcols(object)[["score_pam"]] <- scores
    return(object)
})


#' @rdname addPamScores
#' @export
setMethod("addPamScores", "PairedGuideSet", function(object){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addPamScores(unifiedGuideSet)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})


#' @rdname addPamScores
#' @export
setMethod("addPamScores", "NULL", function(object){
    return(NULL)
})


