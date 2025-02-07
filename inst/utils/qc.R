#' @title Get QC metrics for a functional genomic screen
#' @description Get QC metrics from a SummarizedExperiment object 
#' for each sample: total number of reads, percentage of mapped reads, 
#' mean and median number of reads, number of barcodes with no reads.
#' @param se A SummarizedExperiment object.
#' @param double Should the percentage of mapped reads in the
#' SummarizedExperiment be doubled? \code{TRUE} by default. See details. 
#' @param platform crispr or barcoding. 
#' @return A data.frame containing total number of reads, percentage of
#' reads mapped, median number of reads, mean number of reads, and number of
#' guides with no reads for each sample.
#' @examples 
#' data("seExample")
#' getQCTable(seExample)
#' @export
#' @importFrom matrixStats colMedians
#' @importClassesFrom DelayedArray DelayedMatrix
getQCTable <- function(se, 
                       double = TRUE,
                       platform = c("crispr", "barcoding")
){
    platform <- match.arg(platform)
    if (is(assays(se)[[1]],"DelayedMatrix")) {
        assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    }
    pheno <- colData(se)
    df <- data.frame(samid = colnames(se))
    if ("Sample" %in% colnames(pheno)) {
        df$sample <- pheno$Sample
    }
    df$n.reads <- pheno$nreads
    df$n.mapped <- pheno$nmapped
    df$mappingPerc <- df$n.mapped/df$n.reads * 100
    if (double) {
        df$mappingPerc <- df$mappingPerc * 2
    }
    df$mappingPerc <- round(df$mappingPerc, 1)
    df$meanCount <- round(colMeans(assays(se)[[1]]))
    df$medianCount <- round(colMedians(assays(se)[[1]]))
    df$n.dropouts <- round(colSums(assays(se)[[1]] == 0))

    #QC specific to barcoding:
    if (platform=="barcoding"){
        Y <- as.matrix(assays(se)[[1]])
        df$n.barcodes.cutoff1 <- colSums(Y>0)
        df$n.barcodes.cutoff10 <- colSums(Y>=10)
        df$n.barcodes.cutoff100 <- colSums(Y>=100)
        df$n.dropouts <- NULL
    }
    return(df)
}


#' @title Get QC metrics for a dual guides screen
#' @description Get QC metrics for a dual guides screen
#' @param se A SummarizedExperiment object.
#' @param double Should the percentage of mapped reads in the
#' SummarizedExperiment be doubled? \code{FALSE} by default.
#' 
#' @return A data.frame containing total number of reads, percentage of
#' reads mapped, median number of reads, mean number of reads, and number of
#' guides with no reads for each sample.
#' 
#' @export
getQCTableDual <- function(se, 
                           double=FALSE
){
    if (is(assays(se)[[1]],"DelayedMatrix")) {
        assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    }
    pheno <- colData(se)
    df <- data.frame(samid = colnames(se))
    if ("Sample" %in% colnames(pheno)) {
        df$sample <- pheno$Sample
    }
    cols <- c("none",
              "barcode1.only",
              "barcode2.only", 
              "invalid.pair", 
              "nmapped")
    Z <- pheno[,cols]
    df$n.reads <- rowSums(as.matrix(Z))
    df$n.mapped <- pheno$nmapped
    df$n.invalid.pair <- pheno$invalid.pair
    df$n.barcode1.only <- pheno$barcode1.only
    df$n.barcode2.only <- pheno$barcode2.only
    df$mappingPerc <- df$n.mapped/df$n.reads * 100
    if (double) {
        df$mappingPerc <- df$mappingPerc * 2
    }
    df$mappingPerc <- round(df$mappingPerc, 1)
    df$meanCount   <- round(colMeans(assays(se)[[1]]))
    df$medianCount <- round(colMedians(assays(se)[[1]]))
    df$n.dropouts <- round(colSums(assays(se)[[1]] == 0))
    return(df)
}








#' @title Get QC metrics for cell tracker barcoding experiments
#' @description Get QC metrics from a SummarizedExperiment object 
#' for each sample: total number of reads, percentage of mapped reads, 
#' mean and median number of reads, number of barcodes with no reads.
#' @param object A SummarizedExperiment object
#' @param cutoff Numeric value specifying minimal number of reads required
#'     to call a barcode
#' @importFrom dplyr full_join
#' @importFrom SummarizedExperiment assays colData
#' @export
getQCTableBarcode <- function(object, 
                              cutoff=10){
    qc <- getQCTable(object, platform="barcoding")
    Y  <- assays(object)[[1]]
    ys <- lapply(1:ncol(Y), function(i){
        temp <- Y[,i]
        temp[temp>=cutoff]
    })
    logys <- lapply(ys, function(x) log2(x+1))
    if (!"Sample" %in% colnames(colData(object))){
        names <- colnames(object)
    } else {
        names <- colData(object)$Sample
    }
    
    df    <- data.frame(sample=names, stringsAsFactors=FALSE)
    meds  <- unlist(lapply(logys, median))
    props.cells <- lapply(seq_along(ys), function(i){
        wh <- which(logys[[i]]>= meds[[i]]-2 & logys[[i]]<= meds[[i]]+2)
        a <- sum(2^logys[[i]][wh])
        b <- sum(2^logys[[i]])
        return(a/b*100)
    })
    props.cells <- unlist(props.cells)
    df$perc.reads.4fold <- round(props.cells,1)
    props.cells <- lapply(seq_along(ys), function(i){
        wh <- which(logys[[i]]>= meds[[i]]-4 & logys[[i]]<= meds[[i]]+4)
        a <- sum(2^logys[[i]][wh])
        b <- sum(2^logys[[i]])
        return(a/b*100)
    }) 
    props.cells <- unlist(props.cells)
    df$perc.reads.8fold <- round(props.cells,1)
    qc <- dplyr::full_join(qc,df)
    return(qc)
}


#' @title Get ROC curve for classifying essential genes
#' @description Get receiver operating characteristic (ROC) curve for 
#' classifying essential genes.
#' @param x A numeric vector
#' @param genes A character vector specifying genes associated with x
#' @param species Should be either human or mouse
#' @param rev Should negative values represent signal? TRUE by default. 
#' @param essential.genes.source Source of essential genes
#' @param nonessential.genes.source Source of non-essential genes
#' @param pseudo Should pseudo ROC curves be computed instead?
#'     FALSE by default.
#' @return A data.frame with 2 columns:
#'     \code{specificity} and \code{sensitivity}.
#' @export
#' @importFrom stats approxfun
getEssentialityROC <- function(x,
                               genes,
                               species=c("human", "mouse"),
                               rev=TRUE, 
                               essential.genes.source="hart2014", 
                               nonessential.genes.source="hart2014",
                               pseudo=FALSE
){
    species <- match.arg(species)
    if (rev){
        x <- -x
    }
    if (length(x)!=length(genes)){
        stop("x must have the same length as genes")
    } 
  
    essential_genes <- getEssentials(source=essential.genes.source,
                                     species=species)
    nonessential_genes <- getNonessentials(source=nonessential.genes.source,
                                           species=species)
    isEssential <- genes %in% essential_genes
    # In the case genes are mouse genes:
    if(sum(isEssential)==0){
        stop("[getEssentialityROC] Cannot find essential ",
             "genes among provided genes.")
    }
    if (!pseudo){
        isNonEssential <- genes %in% nonessential_genes
        if (sum(isNonEssential)==0){
            stop("[getEssentialityROC] Cannot find non-essential genes",
                 " among provided genes. Try with pseudo=TRUE.")
        }
    }
    if (pseudo){
        tmp <- .roc_fast(x, as.numeric(isEssential))
    } else {
        wh <- isEssential | isNonEssential
        tmp <- .roc_fast(x[wh], as.numeric(isEssential[wh]))
    }
    out <- data.frame(specificity=tmp$specificities,
                      sensitivity=tmp$sensitivities)
    return(out)
}


.roc_fast <- function(probs, class){
    # Removing non-finite values:
    wh <- which(is.finite(probs))
    if (length(wh)==0){
        stop("[getEssentialityROC] No finite values provided.")
    }
    probs <- probs[wh]
    class <- class[wh]
    # Calculating ROC:
    class_sorted <- class[order(probs, decreasing=T)]
    TPR <- cumsum(class_sorted) / sum(class)
    FPR <- cumsum(class_sorted == 0) / sum(class == 0)
    return(list(sensitivities=TPR,
                specificities=1-FPR))
}




