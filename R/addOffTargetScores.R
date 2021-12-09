#' @title Add CFD and MIT scores to a \linkS4class{GuideSet} object. 
#' @description Add CFD and MIT off-target scores to a \linkS4class{GuideSet}
#'     object. Package \pkg{crisprScore} must be installed.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param max_mm The maximimum number of mismatches between a spacer and
#'     an off-target. Used to select which off-target scores to
#'     include in calculating the aggregated score.
#' @param includeDistance Should a distance penalty for the MIT score be
#'     included? TRUE by default. 
#' @param offset Numeric value specifying an offset to add to the denominator
#'     when calcuting the aggregated score (inverse summation formula). 0 by
#'     default. 
#' 
#' @return \code{guideSet} with off-target score columns appended in
#'     \code{mcols(alignments(guideSet))} and aggregated off-target score
#'     columns appended in \code{mcols(guideSet)}.
#' 
#' @details See \pkg{crisprScore} package for a description of the off-target
#'     scoring methods.
#' 
#' @examples
#' 
#' # Creating a bowtie index:
#' library(Rbowtie)
#' fasta <- system.file(package="crisprDesign", "fasta/chr12.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,
#'                       outdir=outdir,
#'                       force=TRUE,
#'                       prefix="chr12")
#' bowtieIndex <- file.path(outdir, "chr12")
#' 
#' # Adding spacer alignments:
#' guideSetExample <- guideSetExample[1:3]
#' guideSet <- addSpacerAlignments(guideSetExample,
#'                                 bowtie_index=bowtieIndex,
#'                                 n_mismatches=1,
#'                                 txObject=grListExample)
#' 
#' # Adding off-target scores:
#' guideSet <- addOffTargetScores(guideSet)
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{link{addOnTargetScores}} to add on-target scores.
#' 
#' @export
addOffTargetScores <- function(guideSet,
                               max_mm=2,
                               includeDistance=TRUE,
                               offset=0
){
    if (!requireNamespace("crisprScore")){
        message("Please install crisprScore to use 'addOffTargetScores'")
        return(guideSet)
    }
    
    guideSet <- .validateGuideSet(guideSet)
    .checkOffTargetScoresParameters(guideSet=guideSet,
                                    max_mm=max_mm,
                                    includeDistance=includeDistance,
                                    offset=offset)
    
    guideSet <- .addOffTargetScoresToAlignments(
        guideSet=guideSet,
        includeDistance=includeDistance)
    guideSet <- .addOffTargetScoresToGuideSet(guideSet=guideSet,
                                              max_mm=max_mm,
                                              offset=offset)
    return(guideSet)
}



#' @importFrom S4Vectors mcols isTRUEorFALSE
.checkOffTargetScoresParameters <- function(guideSet,
                                            max_mm,
                                            includeDistance,
                                            offset
){
    if (!"alignments" %in% colnames(S4Vectors::mcols(guideSet))){
        stop("Alignments must be added to guideSet prior to calculating ",
             "off-target scores; see ?addSpacerAlignments")
    }
    stopifnot("'max_mm' must be a single non-negative integer value" = {
        is.vector(max_mm, mode="numeric") &&
            length(max_mm) == 1 &&
            max_mm >= 0 &&
            max_mm %% 1 == 0
    })
    stopifnot("'includeDistance' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(includeDistance)
    })
    stopifnot("'offset' must be a single non-negative numeric value" = {
        is.vector(offset, mode="numeric") &&
            length(offset) == 1 &&
            offset >= 0
    })
    invisible(NULL)
}


#' @importFrom utils data
#' @importFrom S4Vectors split mcols<-
.addOffTargetScoresToAlignments <- function(guideSet,
                                            includeDistance=TRUE
){
    crisprNuclease <- crisprNuclease(guideSet)
    utils::data(SpCas9, package="crisprBase", envir=environment())
    stopifnot("Only SpCas9 nuclease is currently supported" = {
        .identicalNucleases(crisprNuclease, SpCas9)
    })
    spacerLengths <- spacerLength(crisprNuclease)
    stopifnot("Only spacers of length 19nt or 20nt are currently supported" = {
        spacerLengths %in% c(19, 20)
    })
    
    alignments <- alignments(guideSet)
    alignmentSpacers  <- as.character(alignments$query)
    alignmentProtospacers <- paste0(as.character(alignments$spacer),
                           as.character(alignments$pam))
    if (spacerLengths == 19){
        alignmentSpacers  <- paste0("G", alignmentSpacers, recycle0=TRUE)
        alignmentProtospacers <- paste0("G", alignmentProtospacers,
                                        recycle0=TRUE)
    }
    score_cfd <- crisprScore::getCFDScores(spacers=alignmentSpacers,
                                           protospacers=alignmentProtospacers)
    alignments$score_cfd <- score_cfd$score
    score_mit <- crisprScore::getMITScores(spacers=alignmentSpacers,
                                           protospacers=alignmentProtospacers,
                                           includeDistance=includeDistance)
    alignments$score_mit <- score_mit$score
    
    guideSetSpacers <- spacers(guideSet, as.character=TRUE)
    alignments <- S4Vectors::split(alignments,
                                   f=factor(alignments$query,
                                            levels=unique(guideSetSpacers)))
    alignments <- alignments[guideSetSpacers]
    names(alignments) <- names(guideSet)
    S4Vectors::mcols(guideSet)[["alignments"]] <- alignments
    
    return(guideSet)
}


#' @importFrom S4Vectors mcols mcols<-
.addOffTargetScoresToGuideSet <- function(guideSet,
                                          max_mm,
                                          offset
){
    alignments <- alignments(guideSet)
    alignments <- as.data.frame(S4Vectors::mcols(alignments),
                                stringsAsFactors=FALSE)
    validMismatchCount <- alignments$n_mismatches <= max_mm
    alignments <- alignments[validMismatchCount, , drop=FALSE]
    alignments <- split(alignments, f=alignments$query)
    .getAggregateScore <- function(score
    ){
        vapply(alignments, function(x){
            x <- x[[score]]
            x <- sum(x, na.rm=TRUE) + offset
            1/x
        }, FUN.VALUE=numeric(1))
    }
    
    cfd <- .getAggregateScore("score_cfd")
    mit <- .getAggregateScore("score_mit")
    
    spacers <- spacers(guideSet, as.character=TRUE)
    cfdIndices <- match(spacers, names(cfd))
    mitIndices <- match(spacers, names(mit))
    S4Vectors::mcols(guideSet)[["score_cfd"]] <- cfd[cfdIndices]
    S4Vectors::mcols(guideSet)[["score_mit"]] <- mit[mitIndices]
    
    return(guideSet)
}
