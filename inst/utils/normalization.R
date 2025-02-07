#' @title Overwrite controls annotation in rowData(se)
#' @description Overwrite controls annotation in rowData(se). This is useful to
#' annotate control features for internal normalization performed by
#' the \code{runVoomScreen} function in \code{gp.sa.screen}.
#' @param se A \code{SummarizedExperiment} object.
#' @param type Should be either \code{NTC} (non-targeting controls), \code{NEG} 
#' (non-essential genes) or \code{Total} (all features are controls). See details. 
#' @param neg.source Source of non-essential genes to be used for NEG normalization.
#' @param species Should be either human or mouse. If \code{NULL}, species will be 
#' inferred from the \code{SummarizedExperiment} metadata.
#' @return A SummarizedExperiment with edited \code{type} column in \code{rowData}
#' @details 
#' * To annotate non-essential genes (NEG), \code{getNonessentials} is used.
#' * To annotate non-targeting controls (NTC), the function assumes that NTC is present 
#' in the column \code{group} in \code{rowData}.
#' * If type==Total, all genes are used for normalization. 
#' @examples 
#' data("seExample")
#' seExample <- overwriteControlType(seExample, type="NEG")
#' @export
overwriteControlType <- function(se,
                                 type=c("NTC", "NEG", "Total"),
                                 neg.source=c("hart2014", "olfactory"),
                                 species=NULL
){
    gene.col <- intersect(c("gene.symbol", "gene_symbol"), colnames(rowData(se)))
    neg.source <- match.arg(neg.source)
    if (is.null(species)){
        species <- .inferSpecies(se)
    }
    neg <- getNonessentials(source=neg.source, species=species)
    #if (species == "mouse") {
    #        neg <- lapply(neg, .simpleCap) %>% unlist
    #}
    type <- match.arg(type)
    rowData(se)$Type <- "Other"
    if (type=="Total"){
        rowData(se)$Type <- "Control"
    }
    if (type=="NTC"){
        if (!"NTC" %in% rowData(se)$group){
            stop("NTCs not found in the `group` field")
        }
        wh <- which(rowData(se)$group=="NTC")
        rowData(se)$Type[wh] <- "Control"
    }
    if (type=="NEG"){
        wh <- which(rowData(se)[[gene.col]] %in% neg)
        rowData(se)$Type[wh] <- "Control"   
    }
    return(se)
}









#' @title Normalization for library size
#' @description Normalization for library size
#' @param se A SummarizedExperiment object.
#' @param fun Median or mean?
#' @return A SummarizedExperiment with counts normalized
#' @examples \dontrun{
#' se <- normalizeTotal(seExample)
#' }
#' @export
#' @importFrom stringr str_extract
normalizeTotal<- function(se,
                          fun=c("median", "mean")
){
    fun <- match.arg(fun)
    factors <- .getNormalizationFactors(se, fun=fun)
    .normalizeWithFactors(se, factors)
}





#' @title Normalization with non-essential genes
#' @description Normalization with non-essential genes.
#' @param se A SummarizedExperiment object.
#' @param gene.col Column in rowData(se) containg gene symbol information.
#' @param species string: 'human' or 'mouse'
#' @param source string: which set of non-essential genes to use?
#' @return A SummarizedExperiment with counts normalized using non-essential genes.
#' @examples \dontrun{
#' se <- normalizeNeg(seExample)
#' }
#' @export
#' @importFrom stringr str_extract
#' 
normalizeNeg <- function(se,
                         gene.col="group",
                         species=NULL, 
                         source=c("hart2014", "olfactory")
){
    source <- match.arg(source)
    if (is.null(species)){
        species <- .inferSpecies(se)
    }
    neg <- getNonessentials(source=source,
                            species=species)
    #if (species == "mouse") {
    #    neg <- lapply(neg, .simpleCap) %>% unlist
    #}
    factors <- .getNormalizationFactors(se,
                                        genes=neg,
                                        gene.col=gene.col)
    .normalizeWithFactors(se, factors)
}




#' @title Normalization with non-targeting controls (NTCs)
#' @description Normalization with non-targeting controls (NTCs).
#' @param se A SummarizedExperiment object.
#' @param gene.col Column in rowData(se) containg NTC information.
#' @return A SummarizedExperiment with counts normalized using NTCs.
#' @examples \dontrun{
#' se <- normalizeNtc(seExample)
#' }
#' @export
normalizeNtc <- function(se,
                         gene.col="group"
){
    factors <- .getNormalizationFactors(se,
                                        genes="NTC",
                                        gene.col=gene.col)
    .normalizeWithFactors(se, factors)
}




#' @importFrom stats median
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importClassesFrom ShadowArray ShadowMatrix
.getNormalizationFactors <- function(se,
                                     genes=NULL,
                                     gene.col="group",
                                     fun=c("median", "mean"), 
                                     na.rm=TRUE
){
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    if (!gene.col %in% colnames(rowData(se))){
        stop("gene.field not found in colnames(rowData(se))")
    }
    fun <- match.arg(fun)
    ann <- rowData(se) 
    #if (!"gene.symbol" %in% colnames(ann)) {
    #    stop("gene.symbol must be in the annotation")
    #}
    if (is.null(genes)) {
        genes <- unique(ann[[gene.col]])
    }
    wh <- which(ann[[gene.col]] %in% genes)
    if (length(wh) != 0) {
        Y <- log2(assays(se)[[1]] + 1)[wh, , drop = FALSE]
    }
    else {
        warning("None of the provided genes for normalization are present in the guide library.")
    }
    ys <- lapply(seq_len(ncol(Y)), function(i) Y[, i])
    ys <- lapply(ys, function(y) y[y != 0])
    if (fun == "median") {
        factors <- vapply(ys, median, na.rm = na.rm, FUN.VALUE=1)
    }
    else {
        factors <- vapply(ys, mean, na.rm = na.rm, FUN.VALUE=1)
    }
    factors <- factors - median(factors, na.rm = na.rm)
    factors
}




#' @importClassesFrom DelayedArray DelayedMatrix
.normalizeWithFactors <- function(se,
                                  factors
){
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    ann <- rowData(se)
    if (ncol(se) != length(factors)) {
        stop("dimensions dont agree")
    }
    Y <- assays(se)[[1]]
    Y <- log2(Y + 1)
    for (i in seq_len(ncol(Y))) {
        temp <- Y[, i]
        temp[temp != 0] <- temp[temp != 0] - factors[i]
        Y[, i] <- temp
    }
    Y <- 2^Y - 1
    Y <- .makePositive(Y)
    assays(se)[[1]] <- Y
    se
}







#' @title Normalization for differental growth in dropout screens
#' @description Normalization for differential growth in dropout screens using temporal
#' quantile normalization. 
#' @param se A SummarizedExperiment object.
#' @param reference.field String specifying the column of \code{colData(se)}
#'     containing the type of each sample (i.e., reference or not).
#' @param reference.level Character vector specifying the reference
#'     levels of the column named by \code{reference.field}.
#' @param replicate.field String specifying the column of \code{colData(se)}
#'     containing the replicate information.
#' @param stratum.field String specifying the column of \code{colData(se)}
#'     containing the variable for stratification
#' @param group.field String specifying the column of \code{colData(se)}
#'     containing the group variable for normalization
#' @param group.reference Character vector specifying the group to
#'     normalize to from the column named by \code{group.field}.
#' @param verbose Should messaged be printed to the console? TRUE by default.
#' @return A SummarizedExperiment with counts normalized for differential growth
#' @import SummarizedExperiment
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @export
normalizeSQN <- function(se,
                         reference.field=NULL,
                         reference.level=NULL,
                         replicate.field=NULL,
                         stratum.field=NULL, #Time
                         group.field=NULL, #Genotype
                         group.reference=NULL,
                         verbose=TRUE
){
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    if (!is.null(reference.field) & !is.null(reference.level)){
        mode <- "ratios"
        if (verbose){
            cat("[normalizeSQN] SQN applied to ratios. \n")
        }
    } else {
        mode <- "counts"
        if (verbose){
            cat("[normalizeSQN] SQN applied to counts. \n")
        }
    }

    pheno <- colData(se)
    strata <- pheno[[stratum.field]]
    stratas <- unique(strata)
    group <- pheno[[group.field]]

    if (mode=="ratios"){
        seRatios <- .getLogRatiosForNormalization(se,
                                                  reference.field=reference.field,
                                                  reference.level=reference.level,
                                                  replicate.field=replicate.field)
        Y <- assays(seRatios)[[1]]
    } else {
        Y <- assays(se)[[1]]
    }

    # Normalizing within stratum:
    for (x in stratas){
        wh <- which(strata==x & group==group.reference)
        if (length(wh)!=0){
            target <- rowMeans(Y[,strata==x & group==group.reference, drop=FALSE], na.rm=TRUE)
        } else {
            # Case where the stratum doesn't have the reference group
            target <- rowMeans(Y[,strata==x, drop=FALSE], na.rm=TRUE)
        }
        wh <- strata==x
        Y[,wh] <- normalize.quantiles.use.target(Y[,wh, drop=FALSE], target=target)
    }

    if (mode=="ratios"){
        assays(seRatios)[[1]] <- Y
        se <- .getInverseLogRatiosForNormalization(se.ratio=seRatios,
                                                   se=se,
                                                   reference.field=reference.field,
                                                   reference.level=reference.level,
                                                   replicate.field=replicate.field)
    } else {
        assays(se)[[1]] <- Y
    }
    return(se)
}




#' @title Normalization for differental growth in dropout screens for focused libraries
#' @description Normalization for differential growth in dropout screens using simple
#' scaling using a set of essential genes
#' @param se A SummarizedExperiment object.
#' @param reference.field String specifying the column of \code{colData(se)} containing the type of each sample (i.e., reference or not).
#' @param reference.level Character vector specifying the reference levels of the column named by \code{reference.field}.
#' @param replicate.field String specifying the column of \code{colData(se)} containing the replicate information.
#' @param stratum.field String specifying the column of \code{colData(se)} containing the variable for stratification
#' @param group.field String specifying the column of \code{colData(se)} containing the group variable for normalization
#' @param group.reference Character vector specifying the group to normalize to from the column named by \code{group.field}.
#' @param gene.field String specifying the field of \code{rowData(se)} that contains the gene identifier for each barcode.
#' @param fun median or mean: function to calculate scaling factors
#' @param species human or mouse
#' @param genes gene names specifying which genes should be used for estimating the scaling factor. If NULL, Hart2014 essential genes will be used. 
#' @return A SummarizedExperiment with counts normalized for differential growth
#' @import SummarizedExperiment
#' @export 
normalizeGrowthRates <- function(se, 
                                 reference.field=NULL,
                                 reference.level=NULL,
                                 replicate.field=NULL,
                                 stratum.field=NULL, #Time
                                 group.field=NULL,
                                 group.reference=NULL,
                                 gene.field=NULL,
                                 fun=c("median", "mean"), 
                                 species=c("human","mouse"),
                                 genes=NULL
){
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    fun <- match.arg(fun)
    species <- match.arg(species)
    if (is.null(genes)){
        genes <- getEssentials(source="hart2014",
                               species=species)
    } 

    seRatios <- .getLogRatiosForNormalization(se,
                                            reference.field=reference.field,
                                            reference.level=reference.level,
                                            replicate.field=replicate.field)
    factors  <- .getDepletionFactors(seRatios,
                                   reference.field=reference.field,
                                   reference.level=reference.level,
                                   stratum.field=stratum.field,
                                   group.field=group.field,
                                   group.reference=group.reference,
                                   gene.field=gene.field,
                                   fun=fun,
                                   genes=genes)
    Y <- assays(seRatios)[[1]]
    Y <- sweep(Y,2,factors, "*")
    assays(seRatios)[[1]] <- Y
    se <- .getInverseLogRatiosForNormalization(se.ratio=seRatios,
                                               se=se,
                                               reference.field=reference.field,
                                               reference.level=reference.level,
                                               replicate.field=replicate.field
                                               )
    return(se)
}




#' @import SummarizedExperiment
#' @importFrom matrixStats colMedians
.getDepletionFactors <- function(se.ratio, 
                                 reference.field=NULL,
                                 reference.level=NULL,
                                 stratum.field=NULL,
                                 group.field=NULL,
                                 group.reference=NULL,
                                 gene.field=NULL,
                                 fun = c("median", "mean"), 
                                 species = c("human","mouse"),
                                 genes=NULL
){
    fun <- match.arg(fun)
    species <- match.arg(species)
    ann   <- rowData(se.ratio)
    pheno <- colData(se.ratio)
    if (!gene.field %in% colnames(ann)) {
        stop("gene.field must be in the annotation")
    }
    if (is.null(genes)){
        genes <- getEssentials(source="hart2014",
                               species=species)
    }
    wh <- which(ann[[gene.field]] %in% genes)
    Y  <- assays(se.ratio)[[1]]
    if (length(wh) != 0) {
        Y <- Y[wh, , drop = FALSE]
    } else {
        warning("None of the provided genes for normalization",
                " are present in the guide library.")
    }
    if (fun == "median") {
        factors <- colMedians(Y, na.rm = TRUE)
    }
    else {
        factors <- colMeans(Y, na.rm = TRUE)
    }

    strata  <- colData(se.ratio)[[stratum.field]]
    stratas <- unique(strata)
    group   <- pheno[[group.field]]
    ref <- pheno[[reference.field]]
    new.factors <- factors
    new.factors[group==group.reference] <- 1
    new.factors[ref==reference.level]   <- 1
    for (i in seq_along(stratas)){
        wh  <- which(strata==stratas[[i]])
        ref    <- which(strata==stratas[[i]] & group==group.reference)
        nonref <- which(strata==stratas[[i]] & group!=group.reference)
        if (length(ref)==0 | length(nonref)==0){
            ratio <- 1
        } else {
            ratio <- abs(median(factors[ref]))/abs(median(factors[nonref]))
        }
        new.factors[nonref] <- ratio
    }
    return(new.factors)
}




#' @import SummarizedExperiment
#' @importClassesFrom DelayedArray DelayedMatrix
.getLogRatiosForNormalization <- function(se,
                                          reference.field=NULL,
                                          reference.level=NULL,
                                          replicate.field=NULL
){

    if (!reference.field %in% colnames(colData(se))){
        stop('reference.field not found in colData(se)')
    }
    group <- colData(se)[[reference.field]]
    if (!reference.level %in% group){
        stop("reference.level not found in reference.field")
    }
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    }  
    Y <- log2(assays(se)[[1]]+1)

    # Replicate specific log-ratios:
    if (!is.null(replicate.field)){
        if (!replicate.field %in% colnames(colData(se))){
           stop('replicate.field not found in colData(se)')
        }
        replicate  <- colData(se)[[replicate.field]]
        replicates <- unique(replicate)

        for (j in seq_along(replicates)){
            temp <- which(replicate==replicates[j])
            ref  <- which(replicate==replicates[j] & group==reference.level)
            Y[, temp] <- Y[,temp] - Y[,ref]
        }
    } else {
          ref  <- which(group==reference.level)
          Y <- Y - rowMeans(Y[,ref, drop=FALSE], na.rm=TRUE)
    }
    assays(se)[[1]] <- Y
    return(se)
}


#' @import SummarizedExperiment
#' @importClassesFrom DelayedArray DelayedMatrix
.getInverseLogRatiosForNormalization <- function(se,
                                                 se.ratio,
                                                 reference.field=NULL,
                                                 reference.level=NULL,
                                                 replicate.field=NULL
){

    if (!reference.field %in% colnames(colData(se))){
        stop('reference.field not found in colData(se)')
    }
    group <- colData(se)[[reference.field]]
    if (!reference.level %in% group){
        stop("reference.level not found in reference.field")
    }
    cond1 <- "DelayedMatrix" %in% class(assays(se)[[1]])
    cond2 <- "ShadowMatrix" %in% class(assays(se)[[1]])
    if(cond1 | cond2){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    Y  <- log2(assays(se)[[1]]+1)
    Yratio <- assays(se.ratio)[[1]]

    if (!is.null(replicate.field)){
        if (!replicate.field %in% colnames(colData(se))){
            stop('replicate.field not found in colData(se)')
        }
        replicate  <- colData(se)[[replicate.field]]
        replicates <- unique(replicate)
      
        for (j in seq_along(replicates)){
            temp <- which(replicate==replicates[j])
            ref  <- which(replicate==replicates[j] & group==reference.level)
            Yratio[, temp] <- Yratio[,temp] + Y[,ref]
        }
    } else {
        ref  <- which(group==reference.level)
        Yratio <- Yratio + rowMeans(Y[,ref, drop=FALSE], na.rm=TRUE)
    }
    Yratio <- 2^Yratio+1
    assays(se.ratio)[[1]] <- Yratio
    return(se.ratio)
}




#' @title Normalization of counts across pooled sublibraries 
#' @description Normalization of counts across pooled sublibraries
#' @param se A SummarizedExperiment object.
#' @param var.lib Variable storing library information. Must be of the same length as the number of rows in `se`
#' @param fun median or mean: function to calculate scaling factors
#' @return A SummarizedExperiment with normalized counts across sublibraries
#' @export
#' @importFrom matrixStats colMedians
normalizeSublibraries <- function(se,
                                  var.lib,
                                  fun=c("median", "mean")
){
    if (length(var.lib)!=nrow(se)){
        stop("var.lib must be of the same length as the number of rows in the SE object.")
    }
    fun <- match.arg(fun)
    Y   <- log2(assays(se)[[1]]+1)
    Y.norm <- Y
    libs <- unique(var.lib)
    factor.matrix <- matrix(0, nrow=length(libs), ncol=ncol(Y))
    for (i in seq_along(libs)){
        lib <- libs[i]
        wh  <- which(var.lib==lib)
        if (fun=="median"){
            factors <- colMedians(Y[wh,], na.rm=TRUE)
        } else {
            factors <- colMeans(Y[wh,], na.rm=TRUE)
        }
        factor.matrix[i, ] <- factors
    }
    medians <- colMedians(factor.matrix)
    factor.matrix <- sweep(factor.matrix, 2, medians, "-")
    for (i in seq_along(libs)){
        lib <- libs[i]
        wh  <- which(var.lib==lib)
        factors <- factor.matrix[i,]
        Y.norm[wh,] <- sweep(Y.norm[wh,],2,factors, "-")
    }
    Y <- 2^Y.norm-1
    if (min(Y)<0){
        Y <- Y+abs(min(Y))
    }
    assays(se)[[1]] <- Y
    return(se)
}

