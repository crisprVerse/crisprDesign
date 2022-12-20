#' @title Add optimal editing site for base editing gRNAs.
#' 
#' @description Add optimal editing site for base editing gRNAs.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param substitution String indicating which substitution
#'     should be used to estimate the optimal editing position. 
#'     E.g. "C2T" will return the optimal editing position for C to T
#'     editing. 
#' @param ... Additional arguments, currently ignored.
#' 
#' @return An updated object with a colum \code{editing_site} added to
#'     \code{mcols(object)}. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @rdname addEditingSites
#' @importFrom crisprBase getEditingSiteFromPamSite
#' @export
setMethod("addEditingSites",
          "GuideSet", 
          function(object,
                   substitution
){
    object <- .validateGuideSet(object)
    nuc <- crisprNuclease(object)
    if (!is(nuc, "BaseEditor")){
        stop("crisprNuclease(object) must be a BaseEditor object.")
    }
    pamSites <- pamSites(object)
    strand <- as.character(strand(object))
    ambiguousStrand <- strand == "*"
    editingSite <- rep(NA, length(object))
    editingSite[!ambiguousStrand] <- getEditingSiteFromPamSite(
        pam_site=pamSites[!ambiguousStrand],
        strand=strand[!ambiguousStrand],
        baseEditor=nuc,
        substitution=substitution)
    mcols(object)$editing_site <- editingSite 
    return(object)
})




#' @rdname addEditingSites
#' @export
setMethod("addEditingSites",
          "PairedGuideSet", 
          function(object,
                   substitution
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addEditingSites(unifiedGuideSet,
                                       substitution=substitution)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addEditingSites
#' @export
setMethod("addEditingSites", "NULL", function(object){
    return(NULL)
})








