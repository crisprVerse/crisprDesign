#' @title Add non-targeting control (NTC) sequences to \linkS4class{GuideSet}
#' @description Add non-targeting control (NTC) sequences to a
#'     \linkS4class{GuideSet} object.
#' 
#' @param object A \linkS4class{GuideSet} or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param ntcs A named character vector of NTC sequences. Sequences must
#'     consist of appropriate DNA or RNA bases, and have the same spacer length
#'     as spacers in \code{object}. Vector names are assigned as IDs and
#'     seqlevels, and must be unique and distinct from IDs and seqnames present
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
#' names(ntcs) <- paste0("ntc_", 1:4)
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
    if (length(ntcs) == 0){
        return(object)
    }
    
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
.validateNtcs <- function(object,
                          ntcs
                          
){
    stopifnot(
        "ntcs must be a named character vector." =
            is.vector(ntcs, mode="character") &&
            !is.null(names(ntcs))
    )
    
    ntcs <- toupper(ntcs)
    
    pattern <- paste0(c("A", "C", "G", "T"), collapse="")
    pattern <- paste0("^[", pattern, "]+$")
    stopifnot(
        "ntcs must exclusively contain DNA bases." =
            all(grepl(pattern, ntcs))
    )
    
    stopifnot(
        "ntcs must have same length as protospacers in object." =
            all(nchar(ntcs) == spacerLength(object))
    )
    
    ntcNames <- names(ntcs)
    stopifnot(
        "ntcs must have unique names." =
            !any(duplicated(ntcNames))
    )
    
    objectIds <- names(object)
    objectSeqlevels <- GenomeInfoDb::seqlevels(object)
    reservedNames <- c(objectIds, objectSeqlevels)
    stopifnot(
        "ntcs must have names distinct from object IDs and seqlevels." =
            length(intersect(ntcNames, reservedNames)) == 0
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
    
    ntcLengths <- rep(spacerLength(object), ntcCount)
    isCircular <- rep(TRUE, ntcCount)
    genome <- rep("ntc", ntcCount)
    
    ntcSeqinfo <- GenomeInfoDb::Seqinfo(seqnames=names(ntcs),
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
    ntc_gs <- GuideSet(
        ids=names(ntcs),
        protospacers=ntcs,
        pams=rep("NA", length(ntcs)),
        seqnames=names(ntcs),
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
    annotationCols <- colnames(S4Vectors::mcols(object))
    baseCols <- colnames(S4Vectors::mcols(ntcGuideSet))
    annotationCols <- setdiff(annotationCols, baseCols)
    for (i in annotationCols){
        dataColumn <- S4Vectors::mcols(object)[[i]]
        if (!is.atomic(dataColumn)){
            S4Vectors::mcols(ntcGuideSet)[[i]] <-
                lapply(seq_along(ntcGuideSet), function(x){
                    x <- S4Vectors::mcols(object)[[i]][0]
                    if (methods::is(x, "GRangesList")){
                        x <- unlist(x)
                    }
                    x
                })
            if (is(dataColumn, "GRangesList")){
                S4Vectors::mcols(ntcGuideSet)[[i]] <-
                    GRangesList(S4Vectors::mcols(ntcGuideSet)[[i]])
            }
        }
    }
    gs <- c(object, ntcGuideSet)
    return(gs)
}
