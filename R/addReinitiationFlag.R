#' @title Add a logical flag for gRNAs leading to potential reinitiation
#' @description Add a logical flag for gRNAs leading to potential reinitiation.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param tx_id String specifiying Ensembl ID for the 
#'     isoform transcript of interested. E.g. "ENST00000311936".
#' @param grnaLocationUpperLimit Integer value specifying the number of 
#'     nucleotides upstream of the start of the CDS in which to search 
#'     for problematic gRNAs. Default value is 150. gRNAS beyond this
#'     value will not be flagged.
#' @param cdsCutoff Numeric value between 0 and 1 to specify the percentage 
#'     of the CDS in which to search for problematic gRNAs.
#'     Default is 0.20. gRNAS beyond this value will not be flagged.
#'     
#' @return The original \code{object} with an appended column 
#'     \code{reinitiationFlag} with logical values. A \code{TRUE} value
#'     indicates a gRNA in proximity of a potential reinitiation site,
#'     and therefore should be avoided. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
addReinitiationFlag <- function(guideSet,
                                tx_id,
                                grnaLocationUpperLimit=150,
                                cdsCutoff=0.20){
    guideSet$reinitiationFlag <- FALSE
    ann <- geneAnnotation(guideSet)
    if (is.null(ann)){
        stop("geneAnnotation must be added first.")
    }
    empty <- sum(ann$tx_id==tx_id, na.rm=TRUE)==0
    if (empty){
        return(guideSet)
    }
    ann <- ann[which(ann$tx_id==tx_id),,drop=FALSE]
    cond1 <- ann$percentCDS<=(cdsCutoff*100)
    cond2 <- ann$aminoAcidIndex < grnaLocationUpperLimit/3
    cond3 <- ann$downtreamATG
    badGuides <- which(cond1 & cond2 & cond3)
    if (length(badGuides)>0){
        badGuides <- rownames(ann)[badGuides]
        wh <- match(badGuides, names(guideSet))
        guideSet$reinitiationFlag[wh] <- TRUE
    }
    return(guideSet)
}

