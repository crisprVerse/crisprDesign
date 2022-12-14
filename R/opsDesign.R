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
#'     either "hamming" (default) or "levenshtein". 
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
                                     dist_method=c("hamming","levenshtein"),
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
                                    dist_method=c("hamming","levenshtein"),
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







#' Design gRNA library for optical pooled screening
#' 
#' Design gRNA library for optical pooled screening
#' 
#' @param df data.frame containing information about candidate
#'     gRNAs from which to build the OPS library. See details. 
#' @param n_guides Integer specifying how many gRNAs per
#'     gene should be selected. 4 by default.
#' @param gene_field String specifying the column in \code{df}
#'     specifying gene names.
#' @param min_dist_edit Integer specifying the minimum distance edit
#'     required for barcodes to be considered dissimilar. Barcodes that
#'     have edit distances less than the min_dist_edit will not be
#'     included in the library. 2 by default. 
#' @param dist_method String specifying distance method. 
#'     Must be either "hamming" (default) or "levenshtein".
#' @param splitByChunks Should distances be calculated in a chunk-wise
#'     manner? FALSE by default. Highly recommended when the set of query
#'     barcodes is large to reduce memory footprint. 
#' 
#' @return A subset of the \code{df} containing the gRNAs
#'     selected for the OPS library. 
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- unique(guideSetExample)
#' guideSet <- addOpsBarcodes(guideSet)
#' guideSet <- guideSet[1:200]
#' 
#' df <- data.frame(ID=names(guideSet),
#'                  spacer=spacers(guideSet, as.character=TRUE),
#'                  opsBarcode=as.character(guideSet$opsBarcode))
#' 
#' # Creating mock gene:
#' df$gene <- rep(paste0("gene",1:10),each=20)
#' df$rank <- rep(1:20,10)
#' opsLib <- designOpsLibrary(df)
#' 
#' @author Jean-Philippe Fortin
#'
#' @export
designOpsLibrary <- function(df,
                             n_guides=4,
                             gene_field="gene",
                             min_dist_edit=2,
                             dist_method=c("hamming","levenshtein"),
                             splitByChunks=FALSE
){
    dist_method <- match.arg(dist_method)
    df <- .validateOpsGrnaInput(df, gene_field)
    

    genes <- unique(df[[gene_field]])
    counts <- rep(0, length(genes))
    grnaList <- list(selected=df[df$rank<=n_guides,,drop=FALSE],
                     candidates=df[df$rank>n_guides,,drop=FALSE],
                     discarded=NULL,
                     genes=genes)
    grnaList <- .initiateOpsLibrary(grnaList,
                                    dist_method=dist_method,
                                    min_dist_edit=min_dist_edit,
                                    splitByChunks=splitByChunks)
    grnaList <- .updateOpsLibrary(grnaList,
                                  gene_field=gene_field,
                                  n_guides=n_guides,
                                  dist_method=dist_method,
                                  min_dist_edit=min_dist_edit,
                                  splitByChunks=splitByChunks)
    out <- .getFinalOpsLibrary(grnaList)
    out <- out[order(out[[gene_field]], out[["rank"]]),,drop=FALSE]
    return(out)
}



#' Update OPS library with additional gRNAs
#' 
#' Update OPS library with additional gRNAs
#' 
#' @param opsLibrary data.frame obtained from \code{designOpsLibrary}.
#' @param df data.frame containing information about additional
#'     candidate gRNAs to add to the OPS library.
#' @param n_guides Integer specifying how many gRNAs per
#'     gene should be selected. 4 by default.
#' @param gene_field String specifying the column in \code{df}
#'     specifying gene names.
#' @param min_dist_edit Integer specifying the minimum distance edit
#'     required for barcodes to be considered dissimilar. Barcodes that
#'     have edit distances less than the min_dist_edit will not be
#'     included in the library. 2 by default. 
#' @param dist_method String specifying distance method. 
#'     Must be either "hamming" (default) or "levenshtein". 
#' @param splitByChunks Should distances be calculated in a chunk-wise
#'     manner? FALSE by default. Highly recommended when the set of query
#'     barcodes is large to reduce memory footprint.
#' 
#' @author Jean-Philippe Fortin
#'
#' @return A data.frame containing the original gRNAs from 
#'    the input \code{opsLibrary} data.frame as well as additional
#'    gRNAs selected from the input data.frame \code{df}.
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- unique(guideSetExample)
#' guideSet <- addOpsBarcodes(guideSet)
#' guideSet1 <- guideSet[1:200]
#' guideSet2 <- guideSet[201:400]
#' 
#' df1 <- data.frame(ID=names(guideSet1),
#'                   spacer=spacers(guideSet1, as.character=TRUE),
#'                   opsBarcode=as.character(guideSet1$opsBarcode))
#' df2 <- data.frame(ID=names(guideSet2),
#'                   spacer=spacers(guideSet2, as.character=TRUE),
#'                   opsBarcode=as.character(guideSet2$opsBarcode))
#' 
#' # Creating mock gene:
#' df1$gene <- rep(paste0("gene",1:10),each=20)
#' df2$gene <- rep(paste0("gene",1:10+10),each=20)
#' df1$rank <- rep(1:20,10)
#' df2$rank <- rep(1:20,10)
#' opsLib <- designOpsLibrary(df1)
#' opsLib <- updateOpsLibrary(opsLib, df2)
#' 
#' @export
updateOpsLibrary <- function(opsLibrary, 
                             df,
                             n_guides=4,
                             gene_field="gene",
                             min_dist_edit=2,
                             dist_method=c("hamming","levenshtein"),
                             splitByChunks=FALSE
){
    dist_method <- match.arg(dist_method)
    df <- .validateOpsGrnaInput(df, gene_field)
    genes <- unique(df[[gene_field]])
    grnaList <- list(selected=opsLibrary,
                     candidates=df,
                     discarded=NULL,
                     genes=genes)
    grnaList <- .updateOpsLibrary(grnaList,
                                  gene_field=gene_field,
                                  n_guides=n_guides,
                                  dist_method=dist_method,
                                  min_dist_edit=min_dist_edit,
                                  splitByChunks=splitByChunks)
    out <- .getFinalOpsLibrary(grnaList)
    out <- out[order(out[[gene_field]], out[["rank"]]),,drop=FALSE]
    return(out)
}






.validateOpsGrnaInput <- function(df,
                                  gene_field
){
    coreCols <- c("ID", "spacer", "opsBarcode")
    if (!all(coreCols %in% colnames(df))){
        diff <- setdiff(coreCols, colnames(df))
        diff <- paste0(diff, collapse=", ")
        stop("The following mandatory columns are missing: ",
             diff, ".")
    }
    if (!gene_field %in% colnames(df)){
        stop("The column specified by gene_field is missing.")
    }
    if (any(duplicated(df[["ID"]]))){
        stop("Some IDs are duplicated.")
    }
    if (!"rank" %in% colnames(df)){
        message("Since 'rank' column is not provided, using ",
                "default order has ranking.")
        df <- .addOpsRank(df, gene_field)
    } else {
        bad <- sum(is.na(df[["rank"]]))!=0
        if (bad){
            stop("Some values are missing in the rank column.")
        }
    }  
    return(df)
}


# Start OPS design with a set of candidate gRNAs
#' @importFrom Matrix rowSums
.initiateOpsLibrary <- function(grnaList,
                                dist_method,
                                min_dist_edit,
                                splitByChunks
){
    selected <- grnaList[["selected"]]
    mat <- getBarcodeDistanceMatrix(queryBarcodes=selected[["opsBarcode"]],
                                    binnarize=TRUE,
                                    dist_method=dist_method,
                                    min_dist_edit=min_dist_edit,
                                    splitByChunks=splitByChunks)
    good <- Matrix::rowSums(mat>0)==0
    # In case all guides are "bad", add first one only:
    if (sum(good)==0){
        good[1] <- TRUE
    }
    grnaList[["selected"]] <- selected[good,]
    grnaList[["candidates"]] <- rbind(grnaList[["candidates"]], 
                                      selected[!good,])
    return(grnaList)
}




# Attemps to add gRNAs to an initial set of chosen gRNAs
# based on the OPS rules
.updateOpsLibrary <- function(grnaList,
                              gene_field,
                              n_guides,
                              dist_method,
                              min_dist_edit,
                              splitByChunks
){
    shouldWeContinue <- TRUE
    while (shouldWeContinue){
        n <- nrow(grnaList[["selected"]])
        grnaList <- .updateOpsLibraryOnce(grnaList,
                                          gene_field=gene_field,
                                          n_guides=n_guides,
                                          dist_method=dist_method,
                                          min_dist_edit=min_dist_edit,
                                          splitByChunks=splitByChunks)
        counts <- table(factor(grnaList[["selected"]][[gene_field]],
                          levels=grnaList[["genes"]]))
        incomplete <- names(which(counts<n_guides))
        remaining  <- grnaList[["candidates"]]
        if (nrow(remaining)==0){
            shouldWeContinue <- FALSE
        } else {
            if (length(incomplete)==0){
                shouldWeContinue <- FALSE
            } else {
                if (sum(incomplete %in% remaining[[gene_field]])==0){
                    shouldWeContinue <- FALSE
                }
            }
        }
    }
    return(grnaList)
}


# Helper function for .updateOpsLibrary
# Attemps to add one gRNA at a time
.updateOpsLibraryOnce <- function(grnaList,
                                  gene_field,
                                  n_guides,
                                  dist_method,
                                  min_dist_edit,
                                  splitByChunks
){
  
    .getCandidates <- function(genes, n){
        cands <- grnaList[["candidates"]]
        if (!"rank" %in% colnames(cands)){
            stop("rank should be a column of the candidate guides.")
        }
        cands <- cands[order(cands$rank),,drop=FALSE]
        cands <- split(cands, f=cands[[gene_field]])
        cands <- cands[genes]
        cands <- lapply(cands, function(x){
            x[seq_len(n),,drop=FALSE]
        })
        cands <- do.call(rbind, cands)
        return(cands)
    }

    .incrementalUpdate <- function(grnaList, cands){
        lib <- grnaList[["selected"]]
        gab <- grnaList[["discarded"]]

        # Improve design by adding first guides that are
        # most divergent:
        dist <- getBarcodeDistanceMatrix(cands[["opsBarcode"]],
                                         dist_method=dist_method,
                                         min_dist_edit=min_dist_edit,
                                         splitByChunks=splitByChunks)
        score <- Matrix::rowSums(dist>0)
        cands <- cands[order(score),,drop=FALSE]

        for (i in seq_len(nrow(cands))){

            barcode <- cands[i,"opsBarcode"]
            dist <- getBarcodeDistanceMatrix(barcode,
                                             lib[["opsBarcode"]],
                                             dist_method=dist_method,
                                             min_dist_edit=min_dist_edit)
            good <- Matrix::rowSums(dist>0)==0
            if (good){
                lib <- rbind(lib, cands[i,,drop=FALSE])
            } else {
                gab <- rbind(gab, cands[i,,drop=FALSE])
            }
        }
        grnaList[["selected"]] <- lib
        grnaList[["discarded"]] <- gab
        return(grnaList)
    }

    .setdiff.grna <- function(set1, set2){
        if (!is.null(set2) & nrow(set2)>0){
            set1 <- set1[!set1[["ID"]] %in% set2[["ID"]],,drop=FALSE]
        }
        return(set1)
    }

    lib <- grnaList[["selected"]]
    geneChoices <- factor(lib[[gene_field]],
                          levels=grnaList[["genes"]])
    geneSets <- split(lib, f=geneChoices)
    ns <- vapply(geneSets, nrow, FUN.VALUE=0)
    incompleteGenes <- names(geneSets)[ns<n_guides]
    cands <- .getCandidates(incompleteGenes, 1)
    if (length(cands)!=0){
        grnaList[["candidates"]] <- .setdiff.grna(grnaList[["candidates"]],
                                                  cands)
        dist <- getBarcodeDistanceMatrix(cands[["opsBarcode"]],
                                         lib[["opsBarcode"]],
                                         dist_method=dist_method,
                                         min_dist_edit=min_dist_edit,
                                         splitByChunks=splitByChunks)
        cands <- cands[Matrix::rowSums(dist)==0,,drop=FALSE]
        grnaList <- .incrementalUpdate(grnaList, cands)
    }
    return(grnaList)
}




.getFinalOpsLibrary <- function(grnaList){
    return(grnaList[["selected"]])
}




#' Validate gRNA library for optical pooled screening
#' 
#' Validate gRNA library for optical pooled screening
#' 
#' @param df data.frame containing information about candidate
#'     gRNAs from which to build the OPS library. See details. 
#' @param min_dist_edit Integer specifying the minimum distance edit
#'     required for barcodes to be considered dissimilar.
#' @param dist_method String specifying distance method. Must be
#'     either "hamming" (default) or "levenshtein". 
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- unique(guideSetExample)
#' guideSet <- addOpsBarcodes(guideSet)
#' df <- data.frame(ID=names(guideSet),
#'                  spacer=spacers(guideSet, as.character=TRUE),
#'                  opsBarcode=as.character(guideSet$opsBarcode))
#' df$gene <- rep(paste0("gene",1:40),each=20)
#' df$rank <- rep(1:20,40)
#' opsLib <- designOpsLibrary(df)
#' opsLib <- validateOpsLibrary(opsLib)
#' 
#' @return The original \code{df} is all checks pass.
#'     Otherwise, a stop error.
#' 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
validateOpsLibrary <- function(df,
                               min_dist_edit=2,
                               dist_method=c("hamming","levenshtein")
){
    dist <- getBarcodeDistanceMatrix(df[["opsBarcode"]],
                                     dist_method=dist_method,
                                     min_dist_edit=min_dist_edit)
    score <- Matrix::rowSums(dist>0)
    if (any(score>0)){
        stop("The library is not valid with the current parameters.")
    }
    return(df)
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


# Get distance costs 
.getCosts <- function(dist_method=c("hamming","levenshtein")
){
    dist_method <- match.arg(dist_method)
    if (dist_method=="hamming"){
        costs  <- list(sub=1, del=1000, ins=1000)
    } else  {
        costs  <- list(sub=1, del=1, ins=1)
    }
    return(costs)
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







