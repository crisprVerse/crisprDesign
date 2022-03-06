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
#'     when calcuting the aggregated score (inverse summation formula).
#'     0 by default.
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
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' fasta <- system.file(package="crisprDesign", "fasta/chr12.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,
#'                       outdir=outdir,
#'                       force=TRUE,
#'                       prefix="chr12")
#' bowtieIndex <- file.path(outdir, "chr12")
#' 
#' # Adding spacer alignments:
#' data(guideSetExample, package="crisprDesign")
#' data(grListExample, package="crisprDesign")
#' guideSetExample <- guideSetExample[1:3]
#' guideSet <- addSpacerAlignments(guideSetExample,
#'                                 aligner="bowtie",
#'                                 aligner_index=bowtieIndex,
#'                                 n_mismatches=1,
#'                                 bsgenome=BSgenome.Hsapiens.UCSC.hg38,
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
    
    guideSet <- .addOffTargetScoresToAlignments(guideSet=guideSet,
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
    utils::data(SpCas9,
                package="crisprBase",
                envir=environment())
    utils::data(CasRx,
                package="crisprBase",
                envir=environment())
    isCas9  <- .identicalNucleases(crisprNuclease,
                                   SpCas9,
                                   checkSpacerLength=FALSE)
    isCasRx <- .identicalNucleases(crisprNuclease,
                                   CasRx,
                                   checkSpacerLength=FALSE)
    if (!isCas9 & !isCasRx){
        stop("Nuclease must be either SpCas9 or CasRx ",
             "for off-target scoring")
    }
    spacerLen <- spacerLength(crisprNuclease)
    if (isCasRx & spacerLen>27){
        stop("For CasRx, spacer length must be at most 27nt ",
             "for off-target scoring.")
    }
    if (isCas9 & spacerLen>20){
        stop("For SpCas9, spacer length must be at most 20nt ",
             "for off-target scoring.")
    }
   
    aln <- alignments(guideSet)
    spacers      <- as.character(aln$spacer)
    protospacers <- as.character(aln$protospacer)
    pams         <- as.character(aln$pam)
    if (isCasRx){
        protospacers <- DNAStringSet(protospacers)
        protospacers <- reverseComplement(protospacers)
        protospacers <- as.character(protospacers)
    }

    if (isCas9 & spacerLen == 19){
        spacers      <- paste0("G", spacers, recycle0=TRUE)
        protospacers <- paste0("G", protospacers, recycle0=TRUE)
    }
    if (isCas9){
        nuclease <- "SpCas9"
    } else if (isCasRx){
        nuclease <- "CasRx"
    }

    if (isCas9 | isCasRx){
        score_cfd <- crisprScore::getCFDScores(spacers=spacers,
                                               protospacers=protospacers,
                                               pams=pams,
                                               nuclease=nuclease)
        aln$score_cfd <- score_cfd$score
    }
    if (isCas9){
        score_mit <- crisprScore::getMITScores(spacers=spacers,
                                               protospacers=protospacers,
                                               pams=pams,
                                               includeDistance=includeDistance)
        aln$score_mit <- score_mit$score
    }

    
    
    guideSetSpacers <- spacers(guideSet, as.character=TRUE)
    aln <- S4Vectors::split(aln,
                            f=factor(aln$spacer,
                                     levels=unique(guideSetSpacers)))
    aln <- aln[guideSetSpacers]
    names(aln) <- names(guideSet)
    S4Vectors::mcols(guideSet)[["alignments"]] <- aln
    
    return(guideSet)
}


#' @importFrom S4Vectors mcols mcols<-
.addOffTargetScoresToGuideSet <- function(guideSet,
                                          max_mm,
                                          offset
){
    aln <- alignments(guideSet)
    aln <- as.data.frame(S4Vectors::mcols(aln),
                         stringsAsFactors=FALSE)
    validMismatchCount <- aln$n_mismatches <= max_mm
    aln <- aln[validMismatchCount, , drop=FALSE]
    aln <- split(aln, f=aln$spacer)

    .getAggregateScore <- function(score){
        vapply(aln, function(x){
            nmm <- x[["n_mismatches"]]
            x <- x[[score]]
            if (sum(nmm==0)==0){
                x <- sum(x, na.rm=TRUE) + 1
            } else {
                x <- sum(x, na.rm=TRUE) + offset
            }
            return(1/x)
        }, FUN.VALUE=numeric(1))
    }
    
    spacers <- spacers(guideSet, as.character=TRUE)
    crisprNuclease <- crisprNuclease(guideSet)
    utils::data(SpCas9,
                package="crisprBase",
                envir=environment())
    utils::data(SpCas9,
                package="crisprBase",
                envir=environment())
    isCas9  <- .identicalNucleases(crisprNuclease,
                                   SpCas9,
                                   checkSpacerLength=FALSE)
    isCasRx <- .identicalNucleases(crisprNuclease,
                                   CasRx,
                                   checkSpacerLength=FALSE)


    #if (isCas9|isCasRx){
    # Not ready to aggregate for CasRx
    # has multiple isoforms create repeat number
    # of on-targets
    if (isCas9){
        cfd <- .getAggregateScore("score_cfd")
        cfdIndices <- match(spacers, names(cfd))
        S4Vectors::mcols(guideSet)[["score_cfd"]] <- cfd[cfdIndices]
    } 
    if (isCas9){
        mit <- .getAggregateScore("score_mit")
        mitIndices <- match(spacers, names(mit))
        S4Vectors::mcols(guideSet)[["score_mit"]] <- mit[mitIndices]
    }
    
    return(guideSet)
}

