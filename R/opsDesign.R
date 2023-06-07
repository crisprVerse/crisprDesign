#' Design gRNA library for optical pooled screening
#' 
#' Design gRNA library for optical pooled screening
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
#' @param n_guides Integer specifying how many gRNAs per
#'     gene should be selected. 4 by default.
#' @param gene_field String specifying the column in \code{mcols(guideSet)}
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
#' @return A subset of the \code{GuideSet} containing the gRNAs
#'     selected for the OPS library. 
#' 
#' @examples
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- unique(guideSetExample)
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
designOpsLibrary <- function(guideSet,
                             n_cycles,
                             rt_direction=c("5prime", "3prime"),
                             n_guides=4,
                             gene_field="gene",
                             min_dist_edit=2,
                             dist_method=c("hamming","levenshtein"),
                             splitByChunks=FALSE
){
    dist_method <- match.arg(dist_method)
    rt_direction <- match.arg(rt_direction)
    gs <- .validateOpsGrnaInput(guideSet, gene_field)
    

    genes <- unique(gs[[gene_field]])
    counts <- rep(0, length(genes))
    grnaList <- list(selected=gs[gs$rank<=n_guides],
                     candidates=gs[gs$rank>n_guides],
                     discarded=NULL,
                     genes=genes)
    grnaList <- .initiateOpsLibrary(grnaList,
                                    n_cycles=n_cycles,
                                    rt_direction=rt_direction,
                                    dist_method=dist_method,
                                    min_dist_edit=min_dist_edit,
                                    splitByChunks=splitByChunks)
    grnaList <- .updateOpsLibrary(grnaList,
                                  n_cycles=n_cycles,
                                  rt_direction=rt_direction,
                                  gene_field=gene_field,
                                  n_guides=n_guides,
                                  dist_method=dist_method,
                                  min_dist_edit=min_dist_edit,
                                  splitByChunks=splitByChunks)
    gs <- .getFinalOpsLibrary(grnaList)
    gs <- gs[order(gs[[gene_field]], gs[["rank"]])]
    return(gs)
}



# Start OPS design with a set of candidate gRNAs
#' @importFrom Matrix rowSums
.initiateOpsLibrary <- function(grnaList,
                                n_cycles,
                                rt_direction=c("5prime", "3prime"),
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
    grnaList[["selected"]] <- selected[good]
    grnaList[["candidates"]] <- c(grnaList[["candidates"]], 
                                  selected[!good])
    return(grnaList)
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
    guideSet <- .validateOpsGrnaInput(guideSet, gene_field)
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












# Attemps to add gRNAs to an initial set of chosen gRNAs
# based on the OPS rules
.updateOpsLibrary <- function(grnaList,
                              n_cycles,
                              rt_direction=c("5prime", "3prime"),
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
                                  n_cycles,
                                  rt_direction=c("5prime", "3prime"),
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
                               n_cycles,
                               rt_direction=c("5prime", "3prime"),
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







