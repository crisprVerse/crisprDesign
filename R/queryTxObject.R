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
    .validateQueryColumn(queryColumn, results)
    .validateQueryValue(queryValue)
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
#' @param tss_window Window size of genomic region flanking TSS to return.
#'     Must be a numeric vector of length 2 consisting of the upstream limit
#'     and downstream limit. Default is \code{c(-500, 500)}, which includes
#'     500bp upstream and downstream of the TSS.
#' 
#' @return A \linkS4class{GRanges} object. Searches yielding no results will
#'     return \code{NULL}.
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @seealso \code{\link{queryTxObject}} for querying gene annotations.
#' 
#' @export
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges promoters
queryTss <- function(tssObject,
                     queryColumn,
                     queryValue,
                     tss_window=NULL # different behavior from documentation
){
    stopifnot("tssObject must be a GRanges object" = {
        is(tssObject, "GRanges")
    })
    .validateQueryColumn(queryColumn, tssObject)
    .validateQueryValue(queryValue)
    hits <- S4Vectors::mcols(tssObject)[[queryColumn]] %in% queryValue
    results <- tssObject[hits]
    if (!is.null(tss_window)){
        tss_window <- .validateTssWindow(tss_window)
        results <- GenomicRanges::promoters(results,
                                            upstream=(-1 * tss_window[1]),
                                            downstream=tss_window[2])
    }
    results <- .nameQueryResults(results)
    return(results)
}







#' @importFrom S4Vectors mcols
.validateQueryColumn <- function(queryColumn,
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



.validateQueryValue <- function(queryValue
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