#' @title To obtain mRNA nucleotide sequences for a transcript
#' 
#' @description  To obtain mRNA nucleotide sequences for a transcript.
#' 
#' @param txids Character vector specifying Ensembl transcript IDs.
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}} for annotating
#'     on-target and off-target alignments using gene annotation.
#' @param bsgenome A \linkS4class{BSgenome} object from which to extract
#'     sequences if \code{x} is a \linkS4class{GRanges} object.
#' 
#' @return DNAStringSet object representing mRNA sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' \dontrun{
#' library(crisprDesignData)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' txObject <- txdb_human
#' txids <- c("ENST00000457313",
#' "ENST00000256078")
#' 
#' out <- getMrnaSequences(txids,
#'     bsgenome=bsgenome,
#'     txObject=txObject)
#' }
#' 
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings DNAStringSet getSeq
#' @importFrom BiocGenerics width
#' @importFrom S4Vectors metadata<-
#' @export
getMrnaSequences <- function(txids,
                             txObject,
                             bsgenome
){
    exons <- txObject[["exons"]]
    exons <- exons[exons$tx_id %in% txids]
    exons <- exons[order(exons$tx_id, exons$exon_rank)]
    if (length(exons)==0){
        stop("txids not found in txObject.")
    }
    mcols(exons)$sequence <- getSeq(bsgenome, exons)


    #Splitting by tx:
    exons <- split(exons, f=exons$tx_id) 
    exons <- as.list(exons)


    sequences <- lapply(seq_along(exons), function(i) {
        sequence <- paste0(exons[[i]]$sequence, collapse='')
        sequence <- DNAStringSet(sequence)
        # Add exon number:
        ns <- width(exons[[i]])
        exon_number <- rep(seq_along(ns), ns)
        mcols(sequence)$exonNumber <- paste0(exon_number, collapse="")
        sequence
    })
    sequences <- do.call(c, sequences)
    names(sequences) <- names(exons)
    return(sequences)
}


#ids <- unique(txObject$exons$tx_id)
#mrnas <- getMrnaSequences(ids,
#                          bsgenome=bsgenome,
#                          txObject=txObject)



#' \dontrun{
#' library(crisprDesignData)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' txObject <- txdb_human
#' txids <- c("ENST00000457313",
#' "ENST00000256078")
#' 
#' out <- getMrnaSequences(txids,
#'     bsgenome=bsgenome,
#'     txObject=txObject)
#' }
#' 
#' 
