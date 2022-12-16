#' @title Add non-targeting control (NTC) sequences to \linkS4class{GuideSet}
#' @description Add non-targeting control (NTC) sequences to a
#'     \linkS4class{GuideSet} object.
#' 
#' @param object A \linkS4class{GuideSet} or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param ntcs A character vector of NTC sequences. Sequences must consist of
#'     appropriate DNA or RNA bases, and have the same spacer length as spacers
#'     in \code{object}.
#' @param ... Additional arguments, currently ignored.
#'     
#' @details NTC sequences are appended as spacers to the \linkS4class{GuideSet}
#'     object. Each NTC sequence is assigned to its own "chromosome" in the
#'     \code{ntc} genome, as reflected in the \linkS4class{Seqinfo} of the
#'     resulting \linkS4class{GuideSet} object. As placeholder values, NTC
#'     ranges are set to \code{0} and strands set to \code{*}.
#'     
#'     All annotation for NTC spacers appended to \code{object} are set to
#'     \code{NA} or empty list elements. To annotate NTC spacers, you must call
#'     the appropriate function after adding NTCs to the \linkS4class{GuideSet}
#'     object.
#' 
#'     
#' @return The original \code{object} with appended \code{ntcs} spacers.
#'     Pre-existing annotation in \code{object} will be set to \code{NA} or
#'     empty list elements for appended NTC spacers.
#' 
#'
#' @examples
#' set.seed(1000)
#' data(guideSetExample, package="crisprDesign")
#' ntcs <- vapply(1:4, function(x){
#'     seq <- sample(c("A", "C", "G", "T"), 20, replace=TRUE)
#'     paste0(seq, collapse="")
#' }, FUN.VALUE=character(1))
#' gs <- addNtcs(guideSetExample, ntcs)
#' gs
#' 
#' 
#' @rdname addNtcs
#' @importFrom GenomeInfoDb seqinfo<-
setMethod("addNtcs", "GuideSet", function(object,
                                          ntcs
){
    object <- .validateGuideSet(object)
    ntcs <- .validateNtcs(object, ntcs)
    
    newSeqinfo <- .updateSeqinfo(object, ntcs)
    ntcGuideSet <- .createNtcGuideSet(object, ntcs, newSeqinfo)
    GenomeInfoDb::seqinfo(object) <- newSeqinfo
    gs <- .mergeNtcGuideSet(object, ntcGuideSet)
    return(gs)
})


#' @rdname addNtcs
#' @export
setMethod("addNtcs",
          "PairedGuideSet", function(object,
                                     ntcs
          ){
              object <- .validateGuideSetOrPairedGuideSet(object)
              unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
              unifiedGuideSet <- addNtcs(unifiedGuideSet,
                                         ntcs=ntcs)
              out <- .addColumnsFromUnifiedGuideSet(object,
                                                    unifiedGuideSet)
              return(out)
})


#' @rdname addNtcs
#' @export
setMethod("addNtcs", "NULL", function(object,
                                      ...){
    return(NULL)
})






#' @importFrom crisprBase targetType
#' @importFrom Biostrings DNA_BASES RNA_BASES
.validateNtcs <- function(object,
                          ntcs
                          
){
    stopifnot(
        "ntcs must be a character vector." =
            is.vector(ntcs, mode="character")
    )
    
    ntcs <- toupper(ntcs)
    
    if (crisprBase::targetType(crisprNuclease(object)) == "DNA"){
        baseSet <- Biostrings::DNA_BASES
    } else {
        baseSet <- Biostrings::RNA_BASES
    }
    pattern <- paste0(baseSet, collapse="")
    pattern <- paste0("^[", pattern, "]+$")
    stopifnot(
        "ntcs must exclusively contain appropriate DNA or RNA bases." =
            all(grepl(pattern, ntcs))
    )
    
    stopifnot(
        "ntcs must have same length as protospacers in object." =
            all(nchar(ntcs) == spacerLength(object))
    )
    
    return(ntcs)
}




#' @importFrom GenomeInfoDb seqlevels seqinfo merge
#' @importClassesFrom GenomeInfoDb Seqinfo
.updateSeqinfo <- function(object,
                           ntcs
){
    ntcCount <- length(ntcs)
    currentSeqlevels <- GenomeInfoDb::seqlevels(object)
    
    ## appended seqlevels must be unique
    ntcNames <- character(0)
    tally <- 0
    while(length(ntcNames) < ntcCount){
        newNames <- paste0("ntc_", seq(from=tally+1, to=tally+ntcCount))
        newNames <- setdiff(newNames, currentSeqlevels)
        ntcNames <- c(ntcNames, newNames)
        tally <- tally + ntcCount
    }
    ntcNames <- ntcNames[seq(ntcCount)] # trim excess names, if any
    
    ntcLengths <- rep(spacerLength(object), ntcCount)
    isCircular <- rep(TRUE, ntcCount)
    genome <- rep("ntc", ntcCount)
    
    ntcSeqinfo <- GenomeInfoDb::Seqinfo(seqnames=ntcNames,
                                        seqlengths=ntcLengths,
                                        isCircular=isCircular,
                                        genome=genome)
    gsSeqinfo <- GenomeInfoDb::seqinfo(object)
    newSeqinfo <- suppressWarnings(
        ## warning due to no sequenve levels in common
        GenomeInfoDb::merge(gsSeqinfo, ntcSeqinfo)
    )
    
    return(newSeqinfo)
}




#' @importFrom GenomeInfoDb seqlevels seqinfo
.createNtcGuideSet <- function(object,
                               ntcs,
                               newSeqinfo
){
    ntcCount <- length(ntcs)
    ntcSeqlevels <- setdiff(GenomeInfoDb::seqlevels(newSeqinfo),
                            GenomeInfoDb::seqlevels(object))
    
    ntc_gs <- GuideSet(
        ids=ntcSeqlevels,
        protospacers=ntcs,
        pams=rep("NA", ntcCount),
        seqnames=ntcSeqlevels,
        pam_site=0,
        strand="*",
        CrisprNuclease=crisprNuclease(object),
        bsgenome=bsgenome(object),
        seqinfo=newSeqinfo
    )
    
    return(ntc_gs)
}




#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
.mergeNtcGuideSet <- function(object,
                              ntcGuideSet
){
    cols <- colnames(S4Vectors::mcols(object))
    for (i in cols){
        data <- S4Vectors::mcols(object)[[i]]
        if (!is.atomic(data) && !methods::is(data, "XStringSet")){
            S4Vectors::mcols(ntcGuideSet)[[i]] <-
                lapply(seq_along(ntcGuideSet), function(x){
                    S4Vectors::mcols(object)[[i]][0]
            })
        }
    }
    gs <- c(object, ntcGuideSet)
    return(gs)
}
