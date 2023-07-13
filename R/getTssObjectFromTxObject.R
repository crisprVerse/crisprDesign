#' @title Extract TSS coordinates from a gene model object
#' 
#' @description Extract TSS coordinates from a gene model object.
#' 
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}} for annotating
#'     on-target and off-target alignments using gene annotation.
#' 
#' @return A GRanges object containing TSS coordinates
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' 
#' data(grListExample, package="crisprDesign")
#' tss <- getTssObjectFromTxObject(grListExample)
#' 
#' @export
getTssObjectFromTxObject <- function(txObject){
    tss <- txObject[["fiveUTRs"]]
    tss <- tss[tss$exon_rank==1]
    cols <- c("tx_id", "gene_id", "gene_symbol")
    mcols(tss) <- mcols(tss)[,cols]
    tss <- tss[!duplicated(mcols(tss)[,cols])]
    tss$promoter <- tss$tx_id
    tss$ID <- paste0(tss$gene_id, "_", tss$promoter, recycle0=TRUE)
    # Making sure we only retain one coordinate:
    if (length(tss) > 0){
        if (as.character(strand(tss))[[1]]=="+"){
            end(tss) <- start(tss)
        } else {
            start(tss) <- end(tss)
        }
    }
    return(tss)
}
