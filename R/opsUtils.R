#' Add optical pooled screening (OPS) barcodes
#' 
#' Add optical pooled screening (OPS) barcodes.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param n_cycles Integer specifying the number of sequencing 
#'     cycles used in the in situ sequencing. This effectively
#'     determines the length of the barcodes to be used for
#'     sequencing.
#' @param rt_direction String specifying from which direction
#'     the reverse transcription of the gRNA spacer sequence
#'     will occur. Must be either "5prime" or "3prime".
#'     "5prime" by default.
#' 
#' @return The original \code{guideSet} object with an additional
#'     column \code{opsBarcode} stored in \code{mcols(guideSet)}.
#'     The column is a \code{DNAStringSet} storing the OPS
#'     barcode. 
#'     
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' guideSetExample <- addOpsBarcodes(guideSetExample)
#' 
#' @author Jean-Philippe Fortin
#' 
#' 
#' @export
addOpsBarcodes <- function(guideSet,
                           n_cycles=9,
                           rt_direction=c("5prime", "3prime")
){
    spacer_len   <- spacerLength(guideSet)
    rt_direction <- match.arg(rt_direction)
    n_cycles <- .validateNCycles(n_cycles, spacer_len)
    spacers <- spacers(guideSet, as.character=TRUE)
    if (rt_direction=="5prime"){
        barcodes <- substr(spacers, 1, n_cycles)
    } else {
        barcodes <- substr(spacers,
                           spacer_len-n_cycles+1,
                           spacer_len)
    }
    barcodes <- DNAStringSet(barcodes)
    mcols(guideSet)[["opsBarcode"]] <- barcodes
    metadata(guideSet)$opsDesign <- list(n_cycles=n_cycles,
                                         rt_direction=rt_direction)
    return(guideSet)
}






#' Extract optical pooled screening (OPS) barcodes
#' 
#' Extract optical pooled screening (OPS) barcodes.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param n_cycles Integer specifying the number of sequencing 
#'     cycles used in the in situ sequencing. This effectively
#'     determines the length of the barcodes to be used for
#'     sequencing.
#' @param rt_direction String specifying from which direction
#'     the reverse transcription of the gRNA spacer sequence
#'     will occur. Must be either "5prime" or "3prime".
#'     "5prime" by default.
#' 
#' @return The original \code{guideSet} object with an additional
#'     column \code{opsBarcode} stored in \code{mcols(guideSet)}.
#'     The column is a \code{DNAStringSet} storing the OPS
#'     barcode. 
#'     
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' barcodes <- extractOpsBarcodes(guideSetExample)
#' 
#' @author Jean-Philippe Fortin
#' 
#' 
#' @export
extractOpsBarcodes <- function(guideSet,
                               n_cycles=9,
                               rt_direction=c("5prime", "3prime")
){
    spacer_len   <- spacerLength(guideSet)
    rt_direction <- match.arg(rt_direction)
    n_cycles <- .validateNCycles(n_cycles, spacer_len)
    spacers <- spacers(guideSet, as.character=TRUE)
    if (rt_direction=="5prime"){
        barcodes <- substr(spacers, 1, n_cycles)
    } else {
        barcodes <- substr(spacers,
                           spacer_len-n_cycles+1,
                           spacer_len)
    }
    return(barcodes)
}





#' Get distance between query and target sets of barcodes
#' 
#' Get distance between query and target sets of barcodes
#' 
#' @param queryBarcodes Character vector of DNA sequences or DNAStringSet. 
#' @param targetBarcodes Optional character vector of DNA sequences 
#'     or DNAStringSet. If NULL, distances will be calculated between
#'     barcodes provided in \code{queryBarcodes}.
#' @param binnarize Should the distance matrix be made binnary?
#'     TRUE by default. See details section. 
#' @param min_dist_edit Integer specifying the minimum distance edit
#'     required for barcodes to be considered dissimilar when
#'     \code{binnarize=TRUE}, ignored otherwise.
#' @param dist_method String specifying distance method. Must be
#'     either "hamming" (default), "levenshtein", or "hybrid".
#'     The hybrid is similar to Levenshtein, but does not allow
#'     for substitutions. 
#' @param ignore_diagonal When \code{targetBarcodes=NULL}, should the
#'     diagonal distances be set to 0 to ignore self distances?
#'     TRUE by default. 
#' @param splitByChunks Should distances be calculated in a chunk-wise
#'     manner? FALSE by default. Highly recommended when the set of query
#'     barcodes is large to reduce memory footprint. 
#' @param n_chunks Integer specifying the number of chunks to be used
#'     when \code{splitByChunks=TRUE}. If NULL (default), number of chunks
#'     will be chosen automatically.
#' 
#' @return A sparse matrix of class \code{dgCMatrix} or \code{dsCMatrix}
#'     in which rows correspond to \code{queryBarcodes} and columns
#'     correspond to \code{targetBarcodes}. If \code{binnarize=TRUE},
#'     a value of 0 indicates that two barcodes have a distance 
#'     greater of equal to \code{min_dist_edit}, otherwise the value 
#'     is 1. If If \code{binnarize=FALSE}, values represent
#'     the actual calculated distances between barcodes.
#' 
#' @examples 
#' data(guideSetExample, package="crisprDesign")
#' guideSetExample <- addOpsBarcodes(guideSetExample)
#' barcodes <- as.character(guideSetExample$opsBarcode)
#' dist <- getBarcodeDistanceMatrix(barcodes, min_dist_edit=2)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
getBarcodeDistanceMatrix <- function(queryBarcodes,
                                     targetBarcodes=NULL,
                                     binnarize=TRUE,
                                     min_dist_edit=NULL,
                                     dist_method=c("hamming",
                                                   "levenshtein",
                                                   "hybrid"),
                                     ignore_diagonal=TRUE,
                                     splitByChunks=FALSE,
                                     n_chunks=NULL
){
    mode <- "notSelf"
    if (is.null(targetBarcodes)){
        targetBarcodes <- queryBarcodes
        mode <- "self"
    }
    if (is(queryBarcodes, "DNAStringSet")){
        queryBarcodes <- as.character(queryBarcodes)
    }
    if (is(targetBarcodes, "DNAStringSet")){
        targetBarcodes <- as.character(targetBarcodes)
    }

    if (is.null(min_dist_edit) & binnarize){
            stop("min_dist_edit must be specified when binnarize=TRUE.")
    }
    if (!splitByChunks | length(queryBarcodes)<=200){
        out <- .getChunkDistanceMatrix(queryBarcodes=queryBarcodes,
                                       targetBarcodes=targetBarcodes,
                                       min_dist_edit=min_dist_edit,
                                       dist_method=dist_method,
                                       binnarize=binnarize)
    } else {
        if (is.null(n_chunks)){
            n_chunks <- ceiling(length(queryBarcodes)/100)
        }
        chunkIndices <- split(seq_along(queryBarcodes),
                              f=as.numeric(cut(seq_along(queryBarcodes),
                                               n_chunks)))
        results <- lapply(seq_len(n_chunks), function(b){
            .getChunkDistanceMatrix(queryBarcodes[chunkIndices[[b]]],
                                    targetBarcodes=targetBarcodes,
                                    min_dist_edit=min_dist_edit,
                                    dist_method=dist_method,
                                    binnarize=binnarize)
                               
        })
        out <- do.call(rbind, results)
    }
    if (mode=="self" & ignore_diagonal){
        diag(out) <- 0
    }
    return(out)
}





# Core function  for getBarcodeDistanceMatrix
#' @importFrom utils adist
#' @importFrom Matrix Matrix
.getChunkDistanceMatrix <- function(queryBarcodes,
                                    targetBarcodes,
                                    min_dist_edit=NULL,
                                    binnarize=TRUE,
                                    dist_method=c("hamming",
                                                  "levenshtein",
                                                  "hybrid"),
                                    sparse=TRUE
){
    dist_method <- match.arg(dist_method)
    costs <- .getCosts(dist_method)
    X <- adist(queryBarcodes,
               targetBarcodes,
               costs=costs)
    if (binnarize){
        if (is.null(min_dist_edit)){
            stop("min_dist_edit must be specified when binnarize=TRUE.")
        }
        X[X<min_dist_edit]  <- 1
        X[X>=min_dist_edit] <- 0
    }
    rownames(X) <- queryBarcodes
    colnames(X) <- targetBarcodes
    if (sparse){
        X <- Matrix(X, sparse=TRUE)
    }
    return(X)
}



.validateOpsGrnaInput <- function(guideSet,
                                  gene_field
){
    df <- mcols(guideSet)
    if (!gene_field %in% colnames(df)){
        stop("The column specified by gene_field is missing.")
    }
    if (any(duplicated(names(df)))){
        stop("Some IDs are duplicated.")
    }
    if (!"rank" %in% colnames(df)){
        stop("'rank' column should be provided")
    } else {
        bad <- sum(is.na(df[["rank"]]))!=0
        if (bad){
            stop("Some values are missing in the rank column.")
        }
    }  
    return(guideSet)
}

# Get distance costs 
.getCosts <- function(dist_method=c("hamming",
                                    "levenshtein",
                                    "hybrid")
){
    dist_method <- match.arg(dist_method)
    if (dist_method=="hamming"){
        costs  <- list(sub=1, del=1000, ins=1000)
    } else if (dist_method=="levenshtein"){
        costs  <- list(sub=1, del=1, ins=1)
    } else {
        costs  <- list(sub=1000, del=1, ins=1)
    }
    return(costs)
}


# Make sure the number of in situ sequencing (ISS) cycles
# is less than the full spacer length
.validateNCycles <- function(n_cycles, spacer_len){
    if (n_cycles>spacer_len){
        stop("n_cycles must be an integer smaller or equal ",
             "to the spacer length.")
    }
    return(n_cycles)
}


.addOpsRank <- function(df, gene_field){
    ids <- df[["ID"]]
    rows <- rownames(df)
    dfs <- split(df, f=df[[gene_field]])
    dfs <- lapply(dfs, function(x){
        x$rank <- seq_len(nrow(x))
        x
    })
    df <- do.call(rbind, dfs)
    df <- df[match(ids, df[["ID"]]),,drop=FALSE]
    rownames(df) <- rows
    return(df)
}





