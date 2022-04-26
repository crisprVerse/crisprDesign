#' @title Convenience function to search for gene coordinates.
#' 
#' @description Convenience function to search for gene coordinates.
#' 
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}}. 
#' @param featureType The genomic feature in \code{txObject} to base your
#'     query on. Must be one of the following: "transcripts", "exons", "cds",
#'     "fiveUTRs", "threeUTRs" or "introns".
#' @param queryColumn Character string specifying the column in 
#'     \code{txObject[[featureType]]} to search for \code{queryValue}(s).
#' @param queryValue Vector specifying the value(s) to search for in
#'     \code{txObject[[featureType]][[queryColumn]]}.
#' 
#' @return A \linkS4class{GRanges} object. Searches yielding no results will
#'     return an empty \linkS4class{GRanges} object.
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @seealso \code{\link{queryTss}} for querying TSS annotations.
#' 
#' @examples
#' 
#' data(grListExample, package="crisprDesign")
#' queryTxObject(grListExample,
#'               featureType="cds",
#'               queryColumn="gene_symbol",
#'               queryValue="IQSEC3")
#' 
#' @export
#' @importFrom S4Vectors mcols
queryTxObject <- function(txObject,
                          featureType=c("transcripts",
                                        "exons",
                                        "cds",
                                        "fiveUTRs",
                                        "threeUTRs",
                                        "introns"),
                          queryColumn,
                          queryValue
){
    txObject <- .validateGRangesList(txObject)
    
    featureType <- match.arg(featureType)
    if (!featureType %in% names(txObject)){  # make check consistent w/ formals
        choices <- paste0(names(txObject), collapse=", ")
        stop("Bad argument for featureType: \"", featureType,
             "\"\n  Choices are: ", choices)
    }
    results <- txObject[[featureType]]
    .checkQueryColumn(queryColumn, results)
    .checkQueryValue(queryValue)
    hits <- S4Vectors::mcols(results)[[queryColumn]] %in% queryValue
    results <- results[hits]
    results <- .nameQueryResults(results)
    return(results)
}








#' @title Convenience function to search for TSS coordinates.
#' 
#' @description Convenience function to search for TSS coordinates.
#' 
#' @param tssObject A \linkS4class{GRanges} containing genomic positions
#'     of transcription starting sites (TSSs).
#' @param queryColumn String specifying which column of \code{mcols(tssObject)}
#'     should be searched for.
#' @param queryValue Character vector specifying the values to search for
#'     in \code{tssObject[[queryColumn]]}.
#' @param tss_window Numeric vector of length 2 establishing the genomic
#'     region to return. The value pair sets the 5 prime and 3 prime limits,
#'     respectively, of the genomic region with respect to the TSS. Use
#'     negative value(s) to set limit(s) upstream of the TSS. Default is
#'     \code{c(-500, 500)}, which includes 500bp upstream and downstream of
#'     the TSS.
#' 
#' @return A \linkS4class{GRanges} object. Searches yielding no results will
#'     return an empty \linkS4class{GRanges} object.
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @examples
#' 
#' data(tssObjectExample, package="crisprDesign")
#' queryTss(tssObjectExample,
#'          queryColumn="gene_symbol",
#'          queryValue="IQSEC3")
#' 
#' 
#' @seealso \code{\link{queryTxObject}} for querying gene annotations.
#' 
#' @export
#' @importFrom methods is
#' @importFrom S4Vectors mcols
queryTss <- function(tssObject,
                     queryColumn,
                     queryValue,
                     tss_window=NULL
){
    stopifnot("tssObject must be a GRanges object" = {
        is(tssObject, "GRanges")
    })
    .checkQueryColumn(queryColumn, tssObject)
    .checkQueryValue(queryValue)
    tss_window <- .validateTssWindow(tss_window)
    hits <- S4Vectors::mcols(tssObject)[[queryColumn]] %in% queryValue
    results <- tssObject[hits]
    results <- .applyTssWindow(results, tss_window)
    results <- .nameQueryResults(results)
    return(results)
}







#' @importFrom S4Vectors mcols
.checkQueryColumn <- function(queryColumn,
                              annotationObject
){
    stopifnot("queryColumn must be a character string" = {
        is.vector(queryColumn, mode="character")
    })
    stopifnot("Only one queryColumn may be used per query" = {
        length(queryColumn) == 1
    })
    if (!queryColumn %in% names(S4Vectors::mcols(annotationObject))){
        stop("queryColumn \"", queryColumn, "\" not found")
    }
    invisible(NULL)
}



.checkQueryValue <- function(queryValue
){
    stopifnot("queryValue must be an atomic vector" = {
        is.null(queryValue) ||
            is.vector(queryValue, mode="logical") ||
            is.vector(queryValue, mode="numeric") ||
            is.vector(queryValue, mode="character")
    })
    invisible(NULL)
}



.nameQueryResults <- function(results
){
    resultsNames <- paste0("region_", seq_along(results))
    names(results) <- resultsNames[seq_along(results)]
    return(results)
}



#' @importFrom BiocGenerics strand
#' @importFrom GenomicRanges promoters
.applyTssWindow <- function(results,
                            tss_window
){
    if (all(tss_window > 0) || all(tss_window < 0)){
        if (all(tss_window > 0)){
            shift <- tss_window[1]
        } else {
            shift <- tss_window[2]
        }
        tss_window <- tss_window - shift
        shift <- rep(shift, length(results))
        minusStrand <- as.character(BiocGenerics::strand(results)) == "-"
        shift[minusStrand] <- -1 * shift[minusStrand]
        results <- shift(results, shift=shift)
    }
    results <- GenomicRanges::promoters(results,
                                        upstream=abs(tss_window[1]),
                                        downstream=abs(tss_window[2]))
    return(results)
}
