#data(tssObjectExample)
#data(grListExample)
#txObject = grListExample
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' 
#' data(grListExample, package="crisprDesign")
#' tss <- getTssObjectFromTxObject(grListExample)
#' 
#' 
#' @export
getTssObjectFromTxObject <- function(txObject){
    tss <- txObject[["fiveUTRs"]]
    tss <- tss[tss$exon_rank==1]
    cols <- c("tx_id", "gene_id", "gene_symbol")
    mcols(tss) <- mcols(tss)[,cols]
    tss <- tss[!duplicated(mcols(tss)[,cols])]
    tss$promoter <- tss$tx_id
    tss$ID <- paste0(tss$gene_symbol, "_", tss$promoter)
    # Making sure we only retain one coordinate:
    if (as.character(strand(tss))[[1]]=="+"){
        end(tss) <- start(tss)
    } else {
        start(tss) <- end(tss)
    }
    return(tss)
}
