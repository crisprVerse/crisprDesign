#' @title Retrieve mRNA sequences
#' 
#' @description A function for retrieving mRNA sequences of select transcripts.
#' 
#' @param txids A character vector of Ensembl transcript IDs. IDs not present
#'     in \code{txObject} will be silently ignored.
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained from \code{\link{TxDb2GRangesList}}. Defines genomic
#'     ranges for \code{txids}.
#' @param bsgenome A \linkS4class{BSgenome} object from which to extract
#'     mRNA sequences.
#' 
#' @return A \linkS4class{DNAStringSet} object of mRNA sequences. Note that
#'     sequences are returned as DNA rather than RNA.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' 
#' if (require("BSgenome.Hsapiens.UCSC.hg38")){
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' data(grListExample)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' txids <- c("ENST00000538872", "ENST00000382841")
#' out <- getMrnaSequences(txids, grListExample, bsgenome)
#' }
#' 
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom BiocGenerics width
#' @importFrom S4Vectors mcols<- split
#' @export
getMrnaSequences <- function(txids,
                             txObject,
                             bsgenome
){
    stopifnot("txids must be a character vector of Ensembl transcript IDs" = {
        is.vector(txids, mode="character")
    })
    txObject <- .validateGRangesList(txObject)
    .isBSGenome(bsgenome)
    
    exons <- txObject[["exons"]]
    exons <- exons[exons$tx_id %in% txids]
    exons <- exons[order(exons$tx_id, exons$exon_rank)]
    S4Vectors::mcols(exons)$sequence <- Biostrings::getSeq(bsgenome, exons)
    exons <- S4Vectors::split(exons, f=exons$tx_id) 
    exons <- as.list(exons)

    sequences <- vapply(exons, function(x){
        paste0(x$sequence, collapse="")
    }, FUN.VALUE=character(1))
    sequences <- Biostrings::DNAStringSet(sequences, use.names=TRUE)
    
    if (length(sequences) > 0){
        exonNumber <- vapply(exons, function(x){
            ns <- BiocGenerics::width(x)
            exonNumber <- rep(seq_along(ns), ns)
            paste0(exonNumber, collapse="")
        }, FUN.VALUE=character(1))
        S4Vectors::mcols(sequences)$exonNumber <- exonNumber
    }
    
    return(sequences)
}
