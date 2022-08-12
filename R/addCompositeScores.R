#' @title Add on-target composite score to a \linkS4class{GuideSet} object.
#' @description Add on-target composite score to a \linkS4class{GuideSet}
#'     object.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param methods Character vector specifying method names for on-target
#'     efficiency prediction algorithms to be used to create the composite
#'     score. Note that the specified scores must be added first to the
#'     object using \code{addOnTargetScores}.
#' @param scoreName String specifying the name of the composite score to
#'     be used as a columm name. Users can choose whatever they like. 
#'     Default is "score_composite".
#' @param ... Additional arguments, currently ignored.
#' 
#' @return \code{guideSet} with column specified by \code{scoreName} 
#'     appended in \code{mcols(guideSet)}.
#' 
#' @details The function creates a composite score across a specified list
#'     of on-target scores by first transforming each individual score into
#'     a rank, and then taking the average rank across all specified methods.
#'     This can improve on-target activity prediction robustness. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addOnTargetScores}} to add on-target scores.
#' 
#' @examples
#' gs <- findSpacers("CCAACATAGTGAAACCACGTCTCTATAAAGAATAAAAAATTAGCCGGGTTA")
#' gs <- addOnTargetScores(gs, methods=c("ruleset1", "crisprater"))
#' gs <- addCompositeScores(gs, methods=c("ruleset1", "crisprater"))
#' 
#' @export
#' @rdname addCompositeScores
setMethod("addCompositeScores", "GuideSet",
    function(object,
             methods=c("azimuth",
                       "ruleset1",
                       "ruleset3",
                       "lindel",
                       "deepcpf1",
                       "deephf",
                       "deepspcas9",
                       "enpamgb",
                       "casrxrf",
                       "crisprater",
                       "crisprscan"),
             scoreName="score_composite"
){
    object <- .validateGuideSet(object)
    crisprNuclease <- crisprNuclease(object)
    methods <- .validateOnTargetScoreMethods(methods=methods,
                                             crisprNuclease=crisprNuclease)
    score_columns <- paste0("score_", methods)
    missing <- methods[!score_columns %in% colnames(mcols(object))]
    if (length(missing)>0){
        stop("The following scores need to be added first with ",
             "addOnTargetScores: ",
             paste0(missing, collapse=","))
    }

    valid <- .validSpacersForOnTargetScores(guideSet=object,
                                            crisprNuclease=crisprNuclease)
    if (!any(valid)){
        warning("No valid chromosome data or PAM sequences; ",
                "cannot calculate on-target scores.",
                immediate.=TRUE)
    }


    score_mat <- mcols(object)[,score_columns][valid,,drop=FALSE]
    score_mat <- apply(score_mat,2, as.numeric)

    .getRankScores <- function(mat){
        ns <- colSums(!is.na(mat))
        maxN <- max(ns)
        factors <- ns/maxN
        ranks <- apply(mat,2, rank, na.last=TRUE, ties.method="first")
        ranks[which(is.na(mat))] <- NA
        ranks <- sweep(ranks,2,factors, "/")
        scores <- rowMeans(ranks, na.rm=TRUE)
        return(scores)
    }
    scores <- .getRankScores(score_mat)
    S4Vectors::mcols(object)[[scoreName]] <- rep(NA,
                                                    length(object))
    S4Vectors::mcols(object)[[scoreName]][valid] <- scores
    return(object)
})




#' @rdname addCompositeScores
#' @export
setMethod("addCompositeScores", "PairedGuideSet",
          function(object,
                   methods=c("azimuth",
                             "ruleset1",
                             "ruleset3",
                             "lindel",
                             "deepcpf1",
                             "deephf",
                             "deepspcas9",
                             "enpamgb",
                             "crisprater",
                             "crisprscan",
                             "casrxrf"),
                   scoreName="score_composite"
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addCompositeScores(unifiedGuideSet,
                                          methods=methods,
                                          scoreName=scoreName)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addCompositeScores
#' @export
setMethod("addCompositeScores", "NULL", function(object){
    return(NULL)
})







