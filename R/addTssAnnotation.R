#' @title Add TSS context annotation to a \linkS4class{GuideSet} object
#' @description Add transcription start site (TSS) context annotation to
#'     spacer sequences stored in a \linkS4class{GuideSet} object.
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param tssObject A \linkS4class{GRanges} object containing TSS coordinates
#'     and annotation.
#' @param anchor A character string specifying which gRNA-specific coordinate
#'     to use (\code{cut_site} or \code{pam_site}) when searching for
#'     overlapping TSS regions.
#' @param tss_window A numeric vector of length 2 establishing the window size
#'     of the genomic region around the TSS to include as the "TSS region".
#'     The values set the upstream and downstream limits, respecitvely. The
#'     default is \code{c(-500, 500)}, which includes 500bp upstream (note the
#'     negative value) and downstream of the TSS.
#' @param ignore.strand If \code{TRUE} (default), includes annotation for
#'     gRNAs irrespective of their target strand. Otherwise, only gRNAs
#'     targeting the gene strand will be annotated.
#' 
#' @return A \linkS4class{GuideSet} object with a \code{tssAnnotation} list
#'     column stored in \code{mcols(guideSet)}. See details section for
#'     descriptions of TSS annotation columns.
#' 
#' @details
#' 
#' \code{mcols(guideSet)[["tssAnnotation"]]} includes all columns from
#'     \code{mcols(tssObject)} in addition to the columns described below.
#' 
#' \itemize{
#' \item \code{chr} — gRNA chromosome name.
#' \item \code{anchor_site} — Genomic coordinate used to search for overlapping
#'     TSS regions.
#' \item \code{strand} — Strand the gRNA is located on.
#' \item \code{tss_id} — The ID for the TSS in \code{tssObject}, if present.
#' \item \code{tss_strand} — Strand the TSS is located on, as provided in
#'     \code{tssObject}
#' \item \code{tss_pos} — Genomic coordinate of the TSS, as provided in
#'     \code{tssObject}.
#' \item \code{dist_to_tss} — Distance (in nucleotides) between the gRNA
#'     \code{anchor_site} and \code{tss_pos}. Negative values indicate
#'     gRNA targets upstream of the TSS.
#' }
#' 
#' @examples 
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- addTssAnnotation(guideSetExample,
#'                              tssObject=tssObjectExample)
#' 
#' # To access TSS annotation:
#' ann <- tssAnnotation(guideSet)
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{\link{addGeneAnnotation}} to add gene annotation, and
#'     \code{\link{tssAnnotation}} to retrieve an existing TSS annotation.
#' 
#' @export
#' @importFrom BiocGenerics rownames
#' @importFrom S4Vectors split mcols<-
addTssAnnotation <- function(guideSet,
                             tssObject,
                             anchor=c("cut_site", "pam_site"),
                             tss_window=NULL,
                             ignore.strand=TRUE
){
    guideSet <- .validateGuideSet(guideSet)
    anchor <- match.arg(anchor)
    tssAnn <- .getTssAnnotation(guideSet=guideSet,
                                tssObject=tssObject,
                                anchor=anchor,
                                tss_window=tss_window,
                                ignore.strand=ignore.strand)
    splitFactor <- factor(BiocGenerics::rownames(tssAnn),
                          levels=names(guideSet))
    tssAnn <- S4Vectors::split(tssAnn, f=splitFactor)
    S4Vectors::mcols(guideSet)[["tssAnnotation"]] <- tssAnn
    return(guideSet)
}



#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics rownames
.getTssAnnotation <- function(guideSet,
                              tssObject,
                              anchor,
                              tss_window,
                              ignore.strand
){
    tssAnn <- .annotateTssOverlaps(guideSet=guideSet,
                                   tssObject=tssObject,
                                   anchor=anchor,
                                   tss_window=tss_window,
                                   ignore.strand=ignore.strand)
    tssAnn <- .addDistToTss(tssAnn)
    tssAnn <- .asDataFrame(tssAnn)
    return(tssAnn)
}



#' @importFrom GenomicRanges GPos promoters findOverlaps
#' @importFrom GenomeInfoDb seqnames seqlevels seqlevels<-
#' @importFrom S4Vectors mcols isTRUEorFALSE queryHits subjectHits
#' @importFrom BiocGenerics strand
.annotateTssOverlaps <- function(guideSet,
                                 tssObject,
                                 anchor,
                                 tss_window,
                                 ignore.strand
                                 
){
    anchor <- .validateAnchor(anchor, guideSet)
    anchorSites <- GenomicRanges::GPos(
        seqnames=GenomeInfoDb::seqnames(guideSet),
        pos=S4Vectors::mcols(guideSet)[[anchor]],
        strand=BiocGenerics::strand(guideSet))
    names(anchorSites)  <- names(guideSet)
    
    tssObject <- .validateTssObject(tssObject)
    tss_window <- .validateTssWindow(tss_window)
    targetAnnotation <- GenomicRanges::promoters(tssObject,
                                                 upstream=(-1*tss_window[1]),
                                                 downstream=tss_window[2])
    
    GenomeInfoDb::seqlevels(anchorSites) <- unique(
        c(GenomeInfoDb::seqlevels(anchorSites),
          GenomeInfoDb::seqlevels(targetAnnotation)))
    stopifnot("'ignore.strand' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(ignore.strand)
    })
    overlaps <- GenomicRanges::findOverlaps(anchorSites,
                                            targetAnnotation,
                                            ignore.strand=ignore.strand)
    tssAnn <- anchorSites[S4Vectors::queryHits(overlaps)]
    tssAnn <- .addTssAnnotationColumns(
        tssAnn=tssAnn,
        targetAnnotation=targetAnnotation,
        targetIndices=S4Vectors::subjectHits(overlaps),
        tssObject=tssObject)
    return(tssAnn)
}



#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics strand colnames start
.addTssAnnotationColumns <- function(tssAnn,
                                     targetAnnotation,
                                     targetIndices,
                                     tssObject
){
    targetAnnCols <- BiocGenerics::colnames(S4Vectors::mcols(targetAnnotation))
    for (i in targetAnnCols){
        targetAnnCol <- S4Vectors::mcols(targetAnnotation)[[i]]
        targetAnnCol <- targetAnnCol[targetIndices]
        if (grepl("^id$", i, ignore.case=TRUE)){
            i <- "tss_id"
        }
        S4Vectors::mcols(tssAnn)[[i]] <- targetAnnCol
    }
    
    targetAnnStrand <- BiocGenerics::strand(targetAnnotation)[targetIndices]
    targetAnnStrand <- as.character(targetAnnStrand)
    S4Vectors::mcols(tssAnn)[["tss_strand"]] <- targetAnnStrand
    
    targetAnnPos <- BiocGenerics::start(tssObject)[targetIndices]
    S4Vectors::mcols(tssAnn)[["tss_pos"]] <- targetAnnPos
    
    return(tssAnn)
}



#' @importFrom GenomicRanges pos
#' @importFrom S4Vectors mcols mcols<-
.addDistToTss <- function(tssAnn
){
    dist_to_tss <- GenomicRanges::pos(tssAnn)-S4Vectors::mcols(tssAnn)$tss_pos
    reverseStrand <- S4Vectors::mcols(tssAnn)$tss_strand == "-"
    dist_to_tss[reverseStrand] <- -1*dist_to_tss[reverseStrand]
    S4Vectors::mcols(tssAnn)[["dist_to_tss"]] <- dist_to_tss
    return(tssAnn)
}
