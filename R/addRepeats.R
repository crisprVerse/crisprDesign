#' @title Annotate a \linkS4class{GuideSet} object with repeat elements
#' @description Add an annotation column to a \linkS4class{GuideSet} object
#'     that identifies spacer sequences overlapping repeat elements.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
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
#' 
#' guideSet <- addRepeats(guideSetExample,
#'                        gr.repeats=grRepeatsExample)
#' 
#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
addRepeats <- function(guideSet,
                       gr.repeats=NULL,
                       ignore.strand=TRUE
){
    guideSet <- .validateGuideSet(guideSet)
    stopifnot("gr.repeats must be a GRanges object" = {
        is(gr.repeats, "GRanges")
    })
    repeatOverlaps <- GenomicRanges::findOverlaps(guideSet,
                                                  gr.repeats,
                                                  ignore.strand=ignore.strand)
    guidesInRepeats <- S4Vectors::queryHits(repeatOverlaps)
    guideSet$inRepeats <- seq_along(guideSet) %in% guidesInRepeats
    return(guideSet)
}


#' @title Remove \linkS4class{GuideSet} gRNAs that overlap repeat elements
#' @description Remove \linkS4class{GuideSet} gRNAs that overlap repeat
#'     elements.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param gr.repeats A \linkS4class{GRanges} object containing
#'     repeat elements regions.
#' @param ignore.strand Should gene strand be ignored when annotating?
#'     TRUE by default. 
#' 
#' @return \code{guideSet} filtered for spacer sequences not overlapping
#'     any repeat elements. An \code{inRepeats} column is also appended in
#'     \code{mcols(guideSet)}.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{link{addRepeats}}.
#' 
#' @examples 
#' 
#' guideSet <- removeRepeats(guideSetExample,
#'                           gr.repeats=grRepeatsExample)
#' 
#' @export
removeRepeats <- function(guideSet,
                          gr.repeats,
                          ignore.strand=TRUE
){
    guideSet <- addRepeats(guideSet,
                           gr.repeats=gr.repeats,
                           ignore.strand=ignore.strand)
    guideSet <- guideSet[!guideSet$inRepeats]
    return(guideSet)
}
