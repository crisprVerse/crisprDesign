#' @title Add a gene-specific exon table to a 
#'      \linkS4class{GuideSet} object.
#' 
#' @title Add a gene-specific exon table to a 
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
#'     output exon table. 
#' @param useConsensusIsoform Should a consensus isoform be used to 
#'     annotate exons? FALSE by default. If TRUE, the isoform constructed
#'     by \code{getConsensusIsoform} will be used. 
#' 
#' @return A \linkS4class{GuideSet} object with a "exonTable" DataFrame
#'     stored in \code{mcols(guideSet)}. The entries in the DataFrame
#'     correspond to the values specified by \code{valueColumn}.
#'     Rows correspond to gRNAs in the GuideSet, columns correspond to 
#'     all exons found in \code{txObject} for gene specified by
#'     \code{gene_id}. 
#' 
#' @examples 
#' if (interactive()){
#'     data(guideSetExample, package="crisprDesign")
#'     data(grListExample, package="crisprDesign")
#'     guideSet <- addGeneAnnotation(guideSetExample,
#'                                   txObject=grListExample)
#'     guideSet <- addExonTable(guideSet,
#'                              gene_id="ENSG00000120645",
#'                              txObject=grListExample)
#' 
#'     guideSet$exonTable
#' }
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addGeneAnnotation}} to add gene annotation and
#'     \code{\link{addTxTable}} to add a transcript table. 
#' 
#' @export
#' @rdname addExonTable
addExonTable <- function(guideSet,
                         gene_id,
                         txObject,
                         valueColumn="percentCDS",
                         useConsensusIsoform=FALSE
){
    if (useConsensusIsoform){
        gs <- addExonTable_consensusIsoform(guideSet=guideSet,
                                            gene_id=gene_id,
                                            txObject=txObject)
    } else {
        gs <- addExonTable_allIsoforms(guideSet=guideSet,
                                       gene_id=gene_id,
                                       txObject=txObject,
                                       valueColumn=valueColumn)
    }
    return(gs)
}



addExonTable_consensusIsoform <- function(guideSet,
                                          gene_id,
                                          txObject
){
    consensus <- getConsensusIsoform(gene_id=gene_id,
                                     txObject=txObject)
    exonids <- unique(consensus$exon_id)
    guides <- names(guideSet)
    out <- matrix(NA,
                  nrow=length(guides),
                  ncol=length(exonids))
    rownames(out) <- guides
    colnames(out) <- exonids
    out <- as.data.frame(out)

    # Looking at overlap:
    gr <- getCutSiteRanges(guideSet,
                           nuclease=crisprNuclease(guideSet))
    df <- findOverlaps(gr, consensus, ignore.strand=TRUE)
    df <- as.data.frame(df)
    for (k in seq_len(nrow(df))){
        out[df[k,1],df[k,2]] <- 1
    }    
    out <- DataFrame(out)
    splitFactor <- factor(BiocGenerics::rownames(out),
                          levels=names(guideSet))
    out <- S4Vectors::split(out, f=splitFactor)
    S4Vectors::mcols(guideSet)[["exonTable"]] <- out
    # mcols(guideSet)$exonTable <- out
    return(guideSet)
}


addExonTable_allIsoforms <- function(guideSet,
                                     gene_id,
                                     txObject,
                                     valueColumn="percentCDS"
){
    stopifnot(.hasGeneAnnotation(guideSet))

    # Let's create an empty data.frame
    exonKey <- .getExon2GeneTable(txObject)
    exonids <- unique(exonKey[exonKey$gene_id==gene_id,,drop=FALSE]$exon_id)
    guides <- names(guideSet)
    out <- matrix(NA, nrow=length(guides), ncol=length(exonids))
    rownames(out) <- guides
    colnames(out) <- exonids
    out <- as.data.frame(out)

    # Now let's get the percentCDS:
    ann <- geneAnnotation(guideSet)
    ann <- ann[ann$gene_id==gene_id,,drop=FALSE]
    if (!valueColumn %in% colnames(ann)){
        stop("valueColumn is not found in geneAnnotation(guideSet).")
    }
    df <- data.frame(ID=rownames(ann),
                     exon_id=ann$exon_id)
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
    S4Vectors::mcols(guideSet)[["exonTable"]] <- out
    # mcols(guideSet)$exonTable <- out
    return(guideSet)
}





