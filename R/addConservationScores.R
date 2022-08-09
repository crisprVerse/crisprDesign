#' @title Add on-target composite score to a \linkS4class{GuideSet} object.
#' 
#' @description Add on-target composite score to a \linkS4class{GuideSet}
#'     object.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param conservationFile String specifing the BigWig file containing
#'     conservation scores.
#' @param nucExtension Number of nucleotides to include on each side of the 
#'     cut site to calculate the conservation score. 2 by default. 
#'     The region will have (2*nucExtension + 1) nucleotides in total.
#' @param fun String specifying the function to use to calculate the final
#'     conservation score in the targeted region. Must be either "max"
#'    (default) or "mean".
#' @param scoreName String specifying the name of the conservation score to
#'     be used as a columm name. Users can choose whatever they like. 
#'     Default is "score_conservation".
#' 
#' @return \code{guideSet} with column specified by \code{scoreName} 
#'     appended in \code{mcols(guideSet)}.
#' 
#' @details The function creates a conservation score for each gRNA
#'     by using the max, or average, conservation score in the genomic 
#'     region where the cut occurs. A BigWig file storing conservation 
#'     stores must be provided. Such files can be downloaded from the
#'     UCSC genome browser. See vignette for more information.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @rdname addConservationScores
#' @importFrom crisprBase getCutSiteRanges
#' @importFrom rtracklayer import.bw
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors mcols<-
#' @importFrom IRanges trim
setMethod("addConservationScores", "GuideSet",
    function(object,
             conservationFile,
             nucExtension=2,
             fun=c("max", "mean"),
             scoreName="score_conservation"
){
    fun <- match.arg(fun)
    object <- .validateGuideSet(object)
    crisprNuclease <- crisprNuclease(object)
    valid <- .validSpacersForOnTargetScores(guideSet=object,
                                            crisprNuclease=crisprNuclease)
    if (!any(valid)){
        warning("No valid chromosome data or PAM sequences; ",
                "cannot calculate conservation scores.",
                immediate.=TRUE)
    }

    if (!file.exists(conservationFile)){
        stop("File specified by 'conservationFile' cannot be found. ")
    }

    # Getting the scores
    gs <- object[valid]
    grCutting <- crisprBase::getCutSiteRanges(gs,
                                              nuclease=crisprNuclease)
    BiocGenerics::start(grCutting) <- BiocGenerics::start(grCutting)-nucExtension
    BiocGenerics::end(grCutting) <- BiocGenerics::end(grCutting)+nucExtension
    grCutting <- IRanges::trim(grCutting)
    scores <- rtracklayer::import.bw(conservationFile,
                                     which=grCutting)
    hits <- IRanges::findOverlaps(grCutting, scores)
    hits <- as.data.frame(hits)
    hits$score <- scores$score[hits$subjectHits]
    scores <- split(hits$score, f=hits$queryHits)
    
    # Aggregating:
    if (fun=="max"){
        scores <- vapply(scores, function(x){
            max(x, na.rm=TRUE)
        }, FUN.VALUE=1)
    } else if (fun=="mean"){
        scores <- vapply(scores, function(x){
            mean(x, na.rm=TRUE)
        }, FUN.VALUE=1)
    }
    scores <- scores[as.character(seq_along(gs))]

    # Adding back to the object:
    S4Vectors::mcols(object)[[scoreName]] <- rep(NA,
                                                 length(object))
    S4Vectors::mcols(object)[[scoreName]][valid] <- scores
    return(object)
})



#library(crisprDesign)
#library(crisprDesignGne)
#library(rtracklayer)
#conservationFile <- getConservationFiles("human")
#guideSet <- guideSetExample
#object=guideSet



#' @rdname addConservationScores
#' @export
setMethod("addConservationScores", "PairedGuideSet",
          function(object,
                   conservationFile,
                   nucExtension=2,
                   fun=c("max", "mean"),
                   scoreName="score_conservation"
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addConservationScores(unifiedGuideSet,
                                             conservationFile=conservationFile,
                                             nucExtension=nucExtension,
                                             fun=fun,
                                             scoreName=scoreName)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addConservationScores
#' @export
setMethod("addConservationScores", "NULL", function(object){
    return(NULL)
})


