#library(GenomicRanges)
#library(crisprDesignData)
#txObject <- txdb_human
#gene_id  <- "ENSG00000049618"

#' @title Get the genomic ranges of a consensus isoform
#' 
#' @description Get the genomic ranges of a consensus isoform. 
#'     The consensus isoform is taken as the union of exons across
#'     all isoforms where overlapping exons are merged to produce
#'     a simplified set through the \code{reduce} method of the
#'     \code{GenomicRanges} package. 
#' 
#' @param gene_id String specifiying Ensembl ID for the 
#'     gene of interest. E.g. "ENSG00000049618". ID must be present in 
#'     \code{txObject$exons$gene_id}. 
#' @param txObject A \linkS4class{TxDb} object or a
#'     \linkS4class{GRangesList} object obtained using
#'     \code{\link{TxDb2GRangesList}} to provide a 
#'     gene model annotation. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @return A \code{GRanges} object. 
#' 
#' @examples
#' data(grListExample)
#' gene_id <- "ENSG00000120645"
#' gr <- getConsensusIsoform(gene_id, grListExample)
#' 
#' @rdname getConsensusIsoform
#' @export
#' 
#' @importFrom BiocGenerics sort
#' @importFrom GenomicRanges reduce
getConsensusIsoform <- function(gene_id,
                                txObject
){
    exons <- txObject$exons
    exons <- exons[exons$gene_id==gene_id]
    if (length(exons)==0){
        stop("Could not find gene ID in txObject$exons$gene_id.")
    }
    exons <- GenomicRanges::reduce(exons)
    exons <- BiocGenerics::sort(exons)
    ss <- as.character(strand(exons))[1]
    if (ss=="+"){
        exons$exon_id <- seq_along(exons)
    } else {
        exons$exon_id <- rev(seq_along(exons))
    }
    exons$exon_id <- paste0(gene_id, "_exon", exons$exon_id)
    names(exons) <- exons$exon_id
    exons$gene_id <- gene_id
    return(exons)
}
