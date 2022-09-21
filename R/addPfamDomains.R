#' @title Add Pfam domains annotation to \linkS4class{GuideSet} object
#' 
#' @description Add Pfam domains annotation to \linkS4class{GuideSet} object.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param pfamTable A \linkS4class{DataFrame} obtained using
#'     \code{\link{preparePfamTable}}.
#' @param ... Additional arguments, currently ignored.
#' 
#' @return An updated object with a colum \code{pfam} added to
#'     \code{geneAnnotation(object)}. 
#' 
#' @details In order to call this function,
#'     the \code{object} must contain a gene annotation
#'     by calling first \code{\link{addGeneAnnotation}}. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso See \code{\link{preparePfamTable}} to prepare the Pfam
#'     domain DataFrame object, and see \code{\link{addGeneAnnotation}}
#'     to add a gene annotation to the \code{object}. 
#' 
#' @rdname addPfamDomains
#' @export
setMethod("addPfamDomains",
          "GuideSet", 
          function(object,
                   pfamTable
){
    object <- .validateGuideSet(object)
    stopifnot(.hasGeneAnnotation(object))
    geneAnn <- geneAnnotation(object)
    bm <- pfamTable[pfamTable$ensembl_transcript_id %in% geneAnn$tx_id, ]

    pfam <- lapply(seq_len(nrow(bm)), function(k){
        txId <- geneAnn$tx_id
        aaIndex <- geneAnn$aminoAcidIndex
        hasMatchingTx <- !is.na(txId) & txId == bm$ensembl_transcript_id[k]
        isNotBeforeDomain <- !is.na(aaIndex) & aaIndex >= bm$pfam_start[k]
        isNotAfterDomain <- !is.na(aaIndex) & aaIndex <= bm$pfam_end[k]
        inPfamDomain <- hasMatchingTx & isNotBeforeDomain & isNotAfterDomain
        spacerPfamDomain <- rep(NA, length(geneAnn))
        spacerPfamDomain[inPfamDomain] <- bm$pfam[k]
        spacerPfamDomain
    })
    pfam <- as.data.frame(pfam, row.names = NULL)
    pfam <- apply(pfam, 1, function(x){
        domains <- x[!is.na(x)]
        domains <- unique(domains)
        if (length(domains) > 0){
            paste0(domains, collapse=';')
        } else {
            NA
        }
    })
    geneAnn[["pfam"]] <- pfam
    splitFactor <- factor(BiocGenerics::rownames(geneAnn),
                          levels=names(object))
    geneAnn <- S4Vectors::split(geneAnn, f=splitFactor)
    S4Vectors::mcols(object)[["geneAnnotation"]] <- geneAnn
    return(object)
})




#' @rdname addPfamDomains
#' @export
setMethod("addPfamDomains",
          "PairedGuideSet", 
          function(object,
                   pfamTable
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addPfamDomains(unifiedGuideSet,
                                      pfamTable=pfamTable)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    
    return(out)
})



#' @rdname addPfamDomains
#' @export
setMethod("addPfamDomains", "NULL", function(object){
    return(NULL)
})








#' @title Obtain Pfam domains from biomaRt
#' 
#' @description Obtain Pfam domains from biomaRt for all
#'    transcripts found in a gene model object.
#' 
#' @param txObject A \linkS4class{TxDb} object or a
#'     \linkS4class{GRangesList} object obtained using
#'     \code{\link{TxDb2GRangesList}} to provide a 
#'     gene model annotation. 
#' @param mart_dataset String specifying dataset to be used by \pkg{biomaRt}
#'     for Pfam domains annotation . E.g. "hsapiens_gene_ensembl".
#' 
#' @return A \linkS4class{DataFrame} object with the following columns:
#' 
#' \itemize{
#' \item \code{ensembl_transcript_id} Ensembl transcript ID.
#' \item \code{pfam} Pfam domain name.
#' \item \code{pfam_start} Start amino acid coordinate of the Pfam domain.
#' \item \code{pfam_end} End amino acid coordinate of the Pfam domain.
#' }
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @examples 
#' data(grListExample, package="crisprDesign")
#' 
#' if (interactive()){
#'     pfamTable <- preparePfamTable(grListExample,
#'                                   mart_dataset="hsapiens_gene_ensembl")
#' }
#' 
#' @export
preparePfamTable <- function(txObject,
                             mart_dataset
){
    if (!requireNamespace("biomaRt", quietly=TRUE)){
        message("Please install the biomaRt package to add Pfam annotation.")
    }
    geneKey <- .getTx2GeneTable(txObject)
    txids <- unique(geneKey$tx_id)
    mart <- biomaRt::useMart('ensembl')
    availableDatasets <- biomaRt::listDatasets(mart)$dataset
    if (is.null(mart_dataset) || !mart_dataset %in% availableDatasets){
        stop("mart_dataset '", mart_dataset, "' is not valid. ",
             "Check available datasets with listDatasets(useMart('ensembl')).")
    }
    mart <- biomaRt::useDataset(mart_dataset, mart=mart)
    attributes <- c('ensembl_transcript_id',
                    'pfam',
                    'pfam_start',
                    'pfam_end')
    filters  <- c('ensembl_transcript_id')
    bm <- biomaRt::getBM(attributes=attributes,
                         filters=filters,
                         values=txids,
                         mart=mart)
    bm <- unique(bm)
    bm <- S4Vectors::DataFrame(bm)
    return(bm)
} 





# library(devtools)
# load_all()
# library(crisprDesignData)
# data(txdb_human)
# txObject=txdb_human
# mart_dataset <- "hsapiens_gene_ensembl"

# pfamTable <- preparePfamTable(txdb_human,
#                               mart_dataset="hsapiens_gene_ensembl")
# data(guideSetExampleFullAnnotation)
# object = guideSetExampleFullAnnotation
# object <- addPfamDomains(object, pfamTable=pfamTable)






