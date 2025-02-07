listGenes <- function(se,
                      gene1.field="gene_symbol_1",
                      gene2.field="gene_symbol_2",
                      sko.position=c("first","second")
){
    sko.position <- match.arg(sko.position)
    if (sko.position=="first"){
        out <- unique(rowData(se)[[gene1.field]])
    } else {
        out <- unique(rowData(se)[[gene2.field]])
    }
    return(out)
}

listPairs <- function(se,
                      gene1.field="gene_symbol_1",
                      gene2.field="gene_symbol_2",
                      exclude.sko=TRUE,
                      sko.control="neg"
){
    gene1 <- rowData(se)[[gene1.field]]
    gene2 <- rowData(se)[[gene2.field]]
    pairs <- unique(paste0(gene1, "_", gene2))
    if (exclude.sko){
        pairs <- pairs[!grepl(sko.control, pairs)]
    }
    return(pairs)
}





#' @title Return all gene pairs from a dual-guides screen
#' 
#' @description Return all gene pairs from a dual-guides screen.
#' 
#' @param se A SummarizedExperiment object.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param type.field String specifying column name in 
#'     \code{ann} containing gRNA pairs class.
#'     "class" by default.
#' @param type.levels Character vector specifying which classes of
#'     gRNA pairs in \code{ann[[type.field]]} should be considered.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
getGenePairs <- function(se, 
                         type.field="class",
                         type.levels=c("target_target",
                                       "or_or",
                                       "eg_eg"),
                         gene1.field="gene_symbol_1",
                         gene2.field="gene_symbol_2"
){
    ann  <- rowData(se)
    temp <- ann[ann[[type.field]] %in% type.levels,]
    temp <- temp[,c(gene1.field, gene2.field)]
    temp <- temp[!duplicated(temp),,drop=FALSE]
    temp <- as.data.frame(temp)
    rownames(temp) <- paste0(temp[,1], "_", temp[,2])
    return(temp)
}







#' @title Subset feature annotation for single-knockout constructs
#' 
#' @description Subset feature annotation for single-knockout constructs
#' 
#' @param se A SummarizedExperiment object.
#' @param sko.position String specifying which gRNA position in the
#'     gRNA pairs should be subsetted. Default is "first".
#' 
#' @return A DataFrame.
#' 
#' @examples
#' ann1 <- getSkoAnnotation(seDualExample)
#' ann2 <- getSkoAnnotation(seDualExample, sko.position="second")
#' 
#' @export
#' @importFrom SummarizedExperiment rowData
getSkoAnnotation <- function(se, 
                             sko.position=c("first", "second")
){
    sko.position <- match.arg(sko.position)
    ann   <- rowData(se)
    cols  <- colnames(ann)
    cols1 <- cols[grepl("_1", cols)]
    cols2 <- cols[grepl("_2", cols)]
    if (sko.position=="first"){
        ann <- ann[ ,cols1, drop=FALSE]
    } else {
        ann <- ann[ ,cols2, drop=FALSE]
    }
    ann <- ann[!duplicated(ann),,drop=FALSE]
    rownames(ann) <- NULL
    colnames(ann) <- gsub("_[1-2]+$", "", colnames(ann))
    return(ann)
}





#' @title Subset a SummarizedExperiment for single-knockout data.
#' 
#' @description Subset a SummarizedExperiment for single-knockout data.
#' 
#' @param se A SummarizedExperiment object.
#' @param assay Numeric value specifying the index of the assay
#'     in \code{assays(se)} to be used. 1 by default.
#' @param sko.position String specifying which gRNA position in the
#'     gRNA pairs should be subsetted. Default is "both".
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param gene.pair.field String specifying colum name in \code{rowData(se)}
#'     containing gene pair name.
#' @param genes An optional character vector of gene names to subset in
#'     \code{gene1.field} or/and \code{gene2.field}.
#' @param aggregate Should data be aggregated at the gene level?
#'     FALSE by default.
#' @param aggregate.fun String specifying function to use for aggregating
#'     data when \code{aggregate=TRUE}. "mean" by default.
#' @param return.matrix Should a matrix be returned instead of a 
#'     SummarizedExperiment? FALSE by default.
#' 
#' @return If \code{return.matrix=FALSE}, a SummarizedExperiment, otherwise
#'     a matrix. 
#' 
#' @examples
#' getSkoData(seDualExample, genes="KRAS")
#' 
#' @export
getSkoData <- function(se,
                       assay=1,
                       sko.position=c("both","first", "second"),
                       sko.control="neg",
                       gene1.field="gene_symbol_1",
                       gene2.field="gene_symbol_2",
                       gene.pair.field="group",
                       genes=NULL,
                       aggregate=FALSE,
                       aggregate.fun=c("mean", "median"),
                       return.matrix=FALSE
){
    aggregate.fun <- match.arg(aggregate.fun)
    sko.position <- match.arg(sko.position)
    if (aggregate & sko.position=="both"){
        stop("aggregate must be false when sko.position='both'.")
    }
    ann <- rowData(se)
    genes1 <- unique(ann[[gene1.field]])
    genes2 <- unique(ann[[gene2.field]])
    genes1 <- setdiff(genes1, sko.control)
    genes2 <- setdiff(genes2, sko.control)
    if (!is.null(genes)){
        genes1 <- intersect(genes1, genes)
        genes2 <- intersect(genes2, genes)
    } else {
        genes <- union(genes1, genes2)
    }
    if (length(genes)==0){
        stop("genes not found for sKO data.")
    }
    pairs1 <- paste0(genes1, "_", sko.control)
    pairs2 <- paste0(sko.control, "_", genes2)
    pairs  <- c(pairs1, pairs2)

    if (sko.position=="first"){
        wh <- which(ann[[gene.pair.field]] %in% pairs1)
    } else if (sko.position=="second"){
        wh <- which(ann[[gene.pair.field]] %in% pairs2)
    } else if (sko.position=="both"){
        wh <- which(ann[[gene.pair.field]] %in% pairs)
    }
    if (length(wh)==0){
        return(NULL)
    } 
    se <- se[wh,]

    if (aggregate){
        se <- .aggregateToGeneLevel(se=se,
                                    sko.position=sko.position,
                                    gene1.field=gene1.field,
                                    gene2.field=gene2.field,
                                    fun=aggregate.fun,
                                    assay=assay)
    }

    if (!return.matrix){
        out <- se
    } else {
        out <- assays(se)[[assay]]
    }
    return(out)
}




.getSkoIndices <- function(se,
                           sko.position=c("both","first", "second"),
                           sko.control="neg",
                           gene1.field="gene_symbol_1",
                           gene2.field="gene_symbol_2",
                           gene.pair.field="group"
){
    sko.position <- match.arg(sko.position)
    ann <- rowData(se)
    genes1 <- unique(ann[[gene1.field]])
    genes2 <- unique(ann[[gene2.field]])
    genes1 <- setdiff(genes1, sko.control)
    genes2 <- setdiff(genes2, sko.control)
    genes  <- union(genes1, genes2)
    pairs1 <- paste0(genes1, "_", sko.control)
    pairs2 <- paste0(sko.control, "_", genes2)
    pairs  <- c(pairs1, pairs2)

    if (sko.position=="first"){
        wh <- which(ann[[gene.pair.field]] %in% pairs1)
    } else if (sko.position=="second"){
        wh <- which(ann[[gene.pair.field]] %in% pairs2)
    } else if (sko.position=="both"){
        wh <- which(ann[[gene.pair.field]] %in% pairs)
    }
    if (length(wh)==0){
        wh <- NULL
    } 
    return(wh)
}



.getDkoIndices <- function(se,
                           sko.control="neg",
                           gene.pair.field="group",
                           gene1.field="gene_symbol_1",
                           gene2.field="gene_symbol_2"
){
    ann <- rowData(se)
    genes1 <- ann[[gene1.field]]
    genes2 <- ann[[gene2.field]]
    targs  <- unique(paste0(genes1, "_", genes2))
    targs <- targs[!grepl(sko.control, targs)]
    wh    <- which(ann[[gene.pair.field]] %in% targs)
    if (length(wh)==0){
        wh <- NULL
    }
    return(wh)
}







#' @title Subset a SummarizedExperiment for double-knockout data.
#' 
#' @description Subset a SummarizedExperiment for double-knockout data.
#' 
#' @param se A SummarizedExperiment object.
#' @param gene.pairs A list of character vectors specifying gene pairs
#'     to be subsetted from \code{gene.pair.field}. Each character vector
#'     must be of length 2 and specifies the 2 gene names of the pair.
#'     Example:  \code{list(c("KRAS", "NRAS"), c("NRAS", "HRAS"))}.
#' @param assay Numeric value specifying the index of the assay
#'     in \code{assays(se)} to be used. 1 by default.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param gene.pair.field String specifying colum name in \code{rowData(se)}
#'     containing gene pair name.
#' @param aggregate Should data be aggregated at the gene level?
#'     FALSE by default.
#' @param aggregate.fun String specifying function to use for aggregating
#'     data when \code{aggregate=TRUE}. "mean" by default.
#' @param return.matrix Should a matrix be returned instead of a 
#'     SummarizedExperiment? FALSE by default.
#' 
#' @return If \code{return.matrix=FALSE}, a SummarizedExperiment, otherwise
#'     a matrix. 
#' 
#' @examples
#' getDkoData(seDualExample,
#'     gene.pairs=list(c("KRAS", "NRAS"), c("KRAS", "HRAS")))
#' 
#' @export 
getDkoData <- function(se,
                       gene.pairs=NULL,
                       assay=1,
                       sko.control="neg",
                       gene1.field="gene_symbol_1",
                       gene2.field="gene_symbol_2",
                       gene.pair.field="group",
                       aggregate=FALSE,
                       aggregate.fun=c("mean", "median"),
                       return.matrix=FALSE
){
    ann <- rowData(se)
    genes1 <- ann[[gene1.field]]
    genes2 <- ann[[gene2.field]]
    choices  <- unique(paste0(genes1, "_", genes2))
    if (!is.null(gene.pairs)){
        targs <- vapply(gene.pairs, paste0, collapse="_", FUN.VALUE="a")
        targs <- intersect(targs, choices)
    } else {
        targs <- choices
    }
    targs <- targs[!grepl(sko.control, targs)]
    wh    <- which(ann[[gene.pair.field]] %in% targs)
    if (length(wh)==0){
        return(NULL)
    }
    se <- se[wh,]

    if (aggregate){
        se <- .aggregateToGenePairLevel(se=se,
                                        gene.pair.field=gene.pair.field,
                                        fun=aggregate.fun,
                                        assay=assay)
    }
    if (!return.matrix){
        out <- se
    } else {
        out <- assays(se)[[assay]]
    }
    return(out)
}




.aggregateToGeneLevel <- function(se,
                                  sko.position=c("first", "second"),
                                  gene1.field="gene_symbol_1",
                                  gene2.field="gene_symbol_2",
                                  fun=c("mean", "median"),
                                  assay=1
){
    fun <- match.arg(fun)
    sko.position <- match.arg(sko.position)
    ann <- rowData(se)
    Y <- assays(se)[[assay]]
    if (sko.position=="first"){
        Ys <- split(as.data.frame(Y),
                    f=ann[[gene1.field]])
    } else {
        Ys <- split(as.data.frame(Y),
                    f=ann[[gene2.field]])
    }
    if (fun=="mean"){
        Ys <- lapply(Ys, colMeans, na.rm=TRUE)
    } else {
        Ys <- lapply(Ys, function(x){
            colMedians(as.matrix(x), na.rm=TRUE)
        })
    }
    Y_new <- do.call(rbind, Ys)
    rownames(Y_new) <- names(Ys)
    colnames(Y_new) <- colnames(se)
    out <- SummarizedExperiment(Y_new)
    return(out)
}


.aggregateToGenePairLevel <- function(se,
                                      gene.pair.field="group",
                                      fun=c("mean", "median"),
                                      assay=1
){
    fun <- match.arg(fun)
    ann <- rowData(se)
    Y <- assays(se)[[assay]]
    Ys <- split(as.data.frame(Y), f=ann[[gene.pair.field]])
    
    if (fun=="mean"){
        Ys <- lapply(Ys, colMeans, na.rm=TRUE)
    } else {
        Ys <- lapply(Ys, function(x){
            colMedians(as.matrix(x), na.rm=TRUE)
        })
    }
    Y_new <- do.call(rbind, Ys)
    rownames(Y_new) <- names(Ys)
    colnames(Y_new) <- colnames(se)
    out <- SummarizedExperiment(Y_new)
    return(out)
}



#' @title Estimate average toxicity of double-knockout constructs.
#' 
#' @description Estimate average toxicity of double-knockout constructs.
#' 
#' @param Y_top Numeric matrix of log2 normalized counts for
#'     later time point samples.
#' @param Y_bottom Numeric matrix of log2 normalized counts for
#'     earlier time point samples.
#' @param ann Feature annotation. Usually extracted using \code{rowData(se)}.
#' @param gene1.field String specifying colum name in \code{ann} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{ann} 
#'     containing gene name for gRNA in position 2.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' Y <- log2(assays(seDualExample)[[1]]+1)
#' pheno <- colData(seDualExample)
#' Y_top <- Y[, pheno$Group=="Day16"]
#' Y_bottom <- Y[, pheno$Group=="Ref"]
#' getDkoToxicity(Y_top,Y_bottom, rowData(seDualExample))
#' 
#' @return A numeric value estimating the average
#'    log-fold change for double-knockout constructs.
#' 
#' @importFrom matrixStats rowMedians
#' @importFrom MASS rlm
#' @export
getDkoToxicity <- function(Y_top,
                           Y_bottom,
                           ann,
                           gene1.field="gene_symbol_1",
                           gene2.field="gene_symbol_2",
                           sko.control="neg"
){
    lfc1 <- rowMedians(Y_bottom, na.rm=TRUE)
    lfc2 <- rowMedians(Y_top, na.rm=TRUE)
    lfc <- lfc2 - lfc1
    names(lfc) <- rownames(Y_top)
    se_lfc <- SummarizedExperiment(as.matrix(lfc),
                                   rowData=ann)
    pairs <- getGenePairs(se_lfc,
                          gene1.field=gene1.field,
                          gene2.field=gene2.field)
    sko1 <- getSkoData(se_lfc,
                       sko.control=sko.control,
                       gene1.field=gene1.field,
                       gene2.field=gene2.field,
                       genes=pairs[[gene1.field]],
                       sko.position="first",
                       return.matrix=TRUE,
                       aggregate=TRUE)
    sko2 <- getSkoData(se_lfc,
                       gene1.field=gene1.field,
                       gene2.field=gene2.field,
                       sko.control=sko.control,
                       genes=pairs[[gene2.field]],
                       sko.position="second",
                       return.matrix=TRUE,
                       aggregate=TRUE)
    dko <- getDkoData(se_lfc,
                      gene1.field=gene1.field,
                      gene2.field=gene2.field,
                      sko.control=sko.control,
                      return.matrix=TRUE,
                      aggregate=TRUE)
    pairs$sko1 <- sko1[match(pairs[[gene1.field]], rownames(sko1))]
    pairs$sko2 <- sko2[match(pairs[[gene2.field]], rownames(sko2))]
    pairs$dko  <- dko[match(rownames(pairs), rownames(dko))]
    dko_expected <- pairs$sko1+pairs$sko2
    model <- rlm(pairs$dko~dko_expected)
    toxicity <- model$coefficients[1]
    names(toxicity) <- NULL
    return(toxicity)
}





.fixAnnotation <- function(se){
    ann <- rowData(se)
    wh <- which(ann$ID_1=="neg")
    ann[wh, "gene_symbol_1"] <- "neg"
    wh <- which(ann$ID_2=="neg")
    ann[wh, "gene_symbol_2"] <- "neg"
    rowData(se) <- ann
    return(se)
}








# #' @export
# extractAllLFCs <- function(se,
#                            gene.pair,
#                            assay=1
# ){
#     x <- strsplit(gene.pair, split = "_")[[1]]
#     gene1 <- x[[1]]
#     gene2 <- x[[2]]
#     sko1 <- .getSkoDataBX(se=se,
#                           gene=gene1,
#                           assay=assay,
#                           sko.position="first")
#     sko2 <- .getSkoDataBX(se=se,
#                           gene=gene2,
#                           assay=assay,
#                           sko.position="second")
#     dko <- .getDkoDataBX(se=se,
#                          gene1=gene1,
#                          gene2=gene2,
#                          assay=assay)
#     y <- c(sko1$logRatios,
#            sko2$logRatios,
#            dko$logRatios)
#     return(y)
# }






# dk.bx   <- .getDkoLogRatiosBX(se=se.new, gene.pairs=gene.pairs)
# dk.bx   <- .getDkoLogRatiosBX(se=se.pilot, gene.pairs=list(c("TP53", "BRIP1")))













# #' @title Calculates log ratios of a dual-guide screen and calculates DLFC between two given genes.
# #' @description generates a list of the top 15 enriched genes for one comparison to a 
# #'        dataframe in the global environment called depleted_hits.
# #' @param se.ratios A Summarized Experiment of logRatios.
# #' @param gene.pair List of gene pairs.
# #' @param sko.control Character string indicating the control (e.g. "neg", "NTC", etc...).
# #' @return dFLC for selected guides
# #' @export
# #' @importFrom matrixStats colMedians
# #' @importFrom SummarizedExperiment colData rowData assays
# #' 
# # this was edited and we should use these 
# getDeltaLFC <- function(se.ratios, 
#                       gene.pair,
#                       sko.control,
#                       ...){
#   ann  <- rowData(se.ratios)
#   sko1 <- getSkoLogRatios(se=se.ratios, genes=gene.pair[[1]][[1]], index=1, sko.control=sko.control, return.matrix=TRUE)
#   sko2 <- getSkoLogRatios(se=se.ratios, genes=gene.pair[[1]][[2]], index=2, sko.control=sko.control, return.matrix=TRUE)
#   dko  <- getDkoLogRatios(se=se.ratios, gene.pairs=gene.pair, return.matrix=TRUE)
#   dflc <- colMedians(dko)-colMedians(sko1)-colMedians(sko2) # take median of the guides for each replicate
#   df   <- data.frame(Group=colData(se.ratios)[,condition.field], dflc=dflc, dko=colMedians(dko), sko1=olMedians(dko))####
#   # df   <- sapply(split(df, f=df$Group), function(x) mean(x$dflc))
#   df   <- vapply(split(df, f=df$Group), function(x) mean(x$dflc), FUN.VALUE=0)
#   df   <- as.data.frame(df)
#   colnames(df) <- "dlfc"
#   df$Group <- rownames(df)
#   df <- df[order(df$Group),]
#   return(df)
# }














