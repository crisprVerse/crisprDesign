#' @title Add a gene-specific transcript table to a 
#'      \linkS4class{GuideSet} object.
#' 
#' @title Add a gene-specific transcript table to a 
#'      \linkS4class{GuideSet} object.
#' 
#' @param guideSet A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param gene_id String specifying gene ID. 
#' @param txObject A \linkS4class{TxDb} object or a
#'     \linkS4class{GRangesList} object obtained using
#'     \code{\link{TxDb2GRangesList}} to provide a 
#'     gene model annotation. 
#' @param valueColumn String specifying column in
#'     \code{geneAnnotation(guideSet)} to use as values in the 
#'     output transcript table. 
#' 
#' @return A \linkS4class{GuideSet} object with a "txTable" DataFrame
#'     stored in \code{mcols(guideSet)}. The entries in the DataFrame
#'     correspond to the values specified by \code{valueColumn}.
#'     Rows correspond to gRNAs in the GuideSet, columns correspond to 
#'     all transcripts found in \code{txObject} for gene specified by
#'     \code{gene_id}.
#' 
#' @examples 
#' if (interactive()){
#'     data(guideSetExample, package="crisprDesign")
#'     data(grListExample, package="crisprDesign")
#'     guideSet <- addGeneAnnotation(guideSetExample,
#'                                   txObject=grListExample)
#'     guideSet <- addTxTable(guideSet,
#'                            gene_id="ENSG00000120645",
#'                            txObject=grListExample)
#' 
#'     guideSet$txTable
#' }
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addGeneAnnotation}} to add gene annotation. 
#' 
#' @rdname addTxTable
#' @importFrom S4Vectors split mcols<-
#' @export
addTxTable <- function(guideSet,
                       gene_id,
                       txObject,
                       valueColumn="percentCDS"
){
    stopifnot(.hasGeneAnnotation(guideSet))

    # Let's create an empty data.frame
    txKey<- .getTx2GeneTable(txObject)
    txids <- txKey[txKey$gene_id==gene_id,,drop=FALSE]$tx_id
    guides <- names(guideSet)
    out <- matrix(NA, nrow=length(guides), ncol=length(txids))
    rownames(out) <- guides
    colnames(out) <- txids
    out <- as.data.frame(out)

    # Now let's get the percentCDS:
    ann <- geneAnnotation(guideSet)
    ann <- ann[ann$gene_id==gene_id,,drop=FALSE]
    if (!valueColumn %in% colnames(ann)){
        stop("valueColumn is not found in geneAnnotation(guideSet).")
    }
    df <- data.frame(ID=rownames(ann),
                     tx_id=ann$tx_id)
    df[[valueColumn]] <- ann[[valueColumn]]
    df <- df[complete.cases(df),,drop=FALSE]

    # Replacing
    for (k in seq_len(nrow(df))){
        out[df[k,1],df[k,2]] <- df[[valueColumn]][k]
    }    
    out <- DataFrame(out)
    splitFactor <- factor(BiocGenerics::rownames(out),
                          levels=names(guideSet))
    out <- S4Vectors::split(out, f=splitFactor)
    S4Vectors::mcols(guideSet)[["txTable"]] <- out
    return(guideSet)
}

