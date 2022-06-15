#' @title Annotate a \linkS4class{GuideSet} object with repeat elements
#' @description Add an annotation column to a \linkS4class{GuideSet} object
#'     that identifies spacer sequences overlapping repeat elements.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param gr.repeats A \linkS4class{GRanges} object containing repeat
#'     elements regions.
#' @param ignore.strand Should gene strand be ignored when annotating?
#'     TRUE by default. 
#' 
#' @return \code{guideSet} with an \code{inRepeats} column appended in
#'     \code{mcols(guideSet)} that signifies whether the spacer sequence
#'     overlaps a repeat element.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{link{removeRepeats}}.
#' 
#' @examples 
#' data(guideSetExample, package="crisprDesign")
#' data(grRepeatsExample, package="crisprDesign")
#' guideSet <- addRepeats(guideSetExample,
#'                        gr.repeats=grRepeatsExample)
#' 
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @rdname addRepeats
setMethod("addRepeats", "GuideSet", function(object,
                                             gr.repeats=NULL,
                                             ignore.strand=TRUE
){
    object <- .validateGuideSet(object)
    stopifnot("gr.repeats must be a GRanges object" = {
        is(gr.repeats, "GRanges")
    })
    repeatOverlaps <- GenomicRanges::findOverlaps(object,
                                                  gr.repeats,
                                                  ignore.strand=ignore.strand)
    guidesInRepeats <- S4Vectors::queryHits(repeatOverlaps)
    object$inRepeats <- seq_along(object) %in% guidesInRepeats
    return(object)
})


#' @export
#' @rdname addRepeats
setMethod("addRepeats", "PairedGuideSet", function(object,
                                                   gr.repeats=NULL,
                                                   ignore.strand=TRUE
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addRepeats(unifiedGuideSet,
                                  gr.repeats=gr.repeats,
                                  ignore.strand=ignore.strand)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addRepeats
#' @export
setMethod("addRepeats", "NULL", function(object){
    return(NULL)
})





#' @title Remove \linkS4class{GuideSet} gRNAs that overlap repeat elements
#' @description Remove \linkS4class{GuideSet} gRNAs that overlap repeat
#'     elements.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param gr.repeats A \linkS4class{GRanges} object containing
#'     repeat elements regions.
#' @param ignore.strand Should gene strand be ignored when annotating?
#'     TRUE by default. 
#' 
#' @return \code{object} filtered for spacer sequences not overlapping
#'     any repeat elements. An \code{inRepeats} column is also appended in
#'     \code{mcols(object)}.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{link{addRepeats}}.
#' 
#' @examples 
#' data(guideSetExample, package="crisprDesign")
#' data(grRepeatsExample, package="crisprDesign")
#' guideSet <- removeRepeats(guideSetExample,
#'                           gr.repeats=grRepeatsExample)
#' 
#' @export
#' @rdname removeRepeats
setMethod("removeRepeats", "GuideSet", function(object,
                                                gr.repeats=NULL,
                                                ignore.strand=TRUE
){
    object <- addRepeats(object,
                         gr.repeats=gr.repeats,
                         ignore.strand=ignore.strand)
    object <- object[!object$inRepeats]
    return(object)
})


#' @export
#' @rdname removeRepeats
setMethod("removeRepeats", "PairedGuideSet", function(object,
                                                      gr.repeats=NULL,
                                                      ignore.strand=TRUE
){
    object <- .validatePairedGuideSet(object)
    gs1 <- addRepeats(first(object),
                      gr.repeats=gr.repeats,
                      ignore.strand=ignore.strand)
    gs2 <- addRepeats(second(object),
                      gr.repeats=gr.repeats,
                      ignore.strand=ignore.strand)
    toKeep <- !gs1$inRepeats & !gs2$inRepeats
    object <- object[toKeep]
    return(object)
})



#' @rdname removeRepeats
#' @export
setMethod("removeRepeats", "NULL", function(object){
    return(NULL)
})

