#' @title Add distance to TSS for a specificed TSS id
#' 
#' @description Add distance to TSS for a specificed TSS id.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param tss_id String specifiying TSS id to calculate the distance.
#'     The column \code{tssAnnotation(object)$tss_id} will be used
#'     to search for the TSS id. 
#' @param ... Additional arguments, currently ignored.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addTssAnnotation}} to add TSS annotation. 
#' 
#' @return A A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object with an additional
#'     metadata column called \code{distance_to_tss} reporting
#'     the distance (in nucleotides) between the TSS position
#'     of the TSS specified by \code{tss_id} and the protospacer
#'     position. The \code{pam_site} coordinate is used as the representative
#'     position of protospacer sequences. 
#' 
#'     Note that a TSS annotation must be available in the \code{object}.
#'     A TSS annotation can be added using \code{addTssAnnotation}. 
#' 
#' @examples
#' data(guideSetExampleFullAnnotation)
#' tss_id <- "ENSG00000120645_P1"
#' gs <- guideSetExampleFullAnnotation
#' gs <- addDistanceToTss(gs, tss_id)
#' 
#' @rdname addDistanceToTss
#' @export
setMethod("addDistanceToTss",
          "GuideSet", 
          function(object,
                   tss_id
){
    if (!"tssAnnotation" %in% colnames(mcols(object))){
        stop("tssAnnotation not found in GuideSet. ",
             "See ?addTssAnnotation to add a TSS annotation first.")
    }
    
    tssAnn <- tssAnnotation(object, unlist=TRUE)
    tssAnn <- tssAnn[tssAnn$tss_id==tss_id,,drop=FALSE]
    wh <- match(names(object), rownames(tssAnn))
    mcols(object)$dist_to_tss <- tssAnn$dist_to_tss[wh]

    return(object)
})


#' @rdname addDistanceToTss
#' @export
setMethod("addDistanceToTss",
          "PairedGuideSet", 
          function(object,
                   tss_id
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addDistanceToTss(unifiedGuideSet,
                                        tss_id=tss_id)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    
    return(out)
})



#' @rdname addIsoformAnnotation
#' @export
setMethod("addDistanceToTss", "NULL", function(object){
    return(NULL)
})



