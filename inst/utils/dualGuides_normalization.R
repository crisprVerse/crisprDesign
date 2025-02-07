#' @title Simple scaling normalization for dual-guides screens.
#' 
#' @description Simple scaling normalization for dual-guides screens.
#' 
#' @param se A SummarizedExperiment object.
#' @param fun String specifying which function should be
#'     used for normalization. "median" by default.
#' @param type.field String specifying column name in 
#'     \code{ann} containing gRNA pairs class.
#'     "class" by default.
#' @param type.levels.dko Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the double
#'     knockout constructs. 
#' @param type.levels.sko1 Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the single
#'     knockout constructs at position 1.
#' @param type.levels.sko2 Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the single
#'     knockout constructs at position 2.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param gene.pair.field String specifying colum name in \code{rowData(se)}
#'     containing gene pair name.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @details If both \code{type.levels.sko1} and \code{type.levels.sko2}
#'     are NULL, normalization for all pairs will be done using 
#'     the pairs specified by \code{type.levels.dko}. If they are not
#'     NULL, normalization will be done separately for single-knockout
#'     (sKO) and double-knockout (dKO) constructs. 
#' 
#' @return A SummarizedExperiment with normalized counts.
#' 
#' @examples
#' normalizeDual(seDualExample)
#' 
#' @importFrom matrixStats colMedians
#' @export
normalizeDual <- function(se,
                          fun=c("median", "mean"),
                          type.field="class",
                          type.levels.dko="or_or",
                          type.levels.sko1=NULL,
                          type.levels.sko2=NULL,
                          gene1.field="gene_symbol_1",
                          gene2.field="gene_symbol_2",
                          gene.pair.field="group",
                          sko.control="neg"
){
    fun <- match.arg(fun)
    if (!type.field %in% colnames(rowData(se))){
        stop("type.field not found in rowData(se)")
    }
    if (is.null(type.levels.sko1) & is.null(type.levels.sko2)){
        cat("[normalizeDual] Normalization will be done jointly for sKO and dKO constructs. \n")
        mode <- "jointly"
    } else {
        cond1 <- is.null(type.levels.sko1) & !is.null(type.levels.sko2)
        cond2 <- !is.null(type.levels.sko1) & is.null(type.levels.sko2)
        if (cond1 | cond2){
            stop("If type.levels.sko1 is provided, type.levels.sko2 must be ",
                 "provided, and vice versa.")
        } else {
            cat("[normalizeDual] Normalization will be done separately for sKO and dKO constructs. \n")
            mode <- "separately"
        }
    }

    # Extracting data:
    assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    Y <- assays(se)[[1]]
    Y <- log2(Y+1)

    if (mode=="jointly"){
        wh <- which(rowData(se)[[type.field]]==type.levels.dko)
        if (length(wh)==0){
            stop("None of the type.levels.dko are found in the specified type.field column.")
        }
        if (fun=="median"){
            factors <- colMedians(as.matrix(Y)[wh,,drop=FALSE], na.rm=TRUE)
        } else if (fun=="mean"){
            factors <- colMeans(as.matrix(Y)[wh,,drop=FALSE], na.rm=TRUE)
        }
        factors <- factors-median(factors)
        Y <- sweep(Y,2,factors, "-")
    } else {
        wh.dko <- which(rowData(se)[[type.field]]==type.levels.dko)
        wh.sko1 <- which(rowData(se)[[type.field]]==type.levels.sko1)
        wh.sko2 <- which(rowData(se)[[type.field]]==type.levels.sko2)
        if (length(wh.dko)==0){
            stop("None of the type.levels.dko are found in the specified type.field column.")
        }
        if (length(wh.sko1)==0){
            stop("None of the type.levels.sko1 are found in the specified type.field column.")
        }
        if (length(wh.sko2)==0){
            stop("None of the type.levels.sko2 are found in the specified type.field column.")
        }
        if (fun=="median"){
            factors.dko  <- colMedians(as.matrix(Y)[wh.dko,,drop=FALSE], na.rm=TRUE)
            factors.sko1 <- colMedians(as.matrix(Y)[wh.sko1,,drop=FALSE], na.rm=TRUE)
            factors.sko2 <- colMedians(as.matrix(Y)[wh.sko2,,drop=FALSE], na.rm=TRUE)
        } else if (fun=="mean"){
            factors.dko  <- colMeans(as.matrix(Y)[wh.dko,,drop=FALSE], na.rm=TRUE)
            factors.sko1 <- colMeans(as.matrix(Y)[wh.sko1,,drop=FALSE], na.rm=TRUE)
            factors.sko2 <- colMeans(as.matrix(Y)[wh.sko2,,drop=FALSE], na.rm=TRUE)
        }
        factors.dko <- factors.dko-median(factors.dko)
        factors.sko1 <- factors.sko1-median(factors.sko1)
        factors.sko2 <- factors.sko2-median(factors.sko2)
        sko1.indices <- .getSkoIndices(se=se,
                                       sko.position="first",
                                       sko.control=sko.control,
                                       gene1.field=gene1.field,
                                       gene2.field=gene2.field,
                                       gene.pair.field=gene.pair.field)
        sko2.indices <- .getSkoIndices(se=se,
                                       sko.position="second",
                                       sko.control=sko.control,
                                       gene1.field=gene1.field,
                                       gene2.field=gene2.field,
                                       gene.pair.field=gene.pair.field)
        dko.indices <- .getDkoIndices(se=se,
                                      sko.control=sko.control,
                                      gene1.field=gene1.field,
                                      gene2.field=gene2.field,
                                      gene.pair.field=gene.pair.field)
        Y[sko1.indices,] <- sweep(Y[sko1.indices,,drop=FALSE],2,factors.sko1, "-")
        Y[sko2.indices,] <- sweep(Y[sko2.indices,,drop=FALSE],2,factors.sko2, "-")
        Y[dko.indices,]  <- sweep(Y[dko.indices,,drop=FALSE],2,factors.dko, "-")
    }

    # Transforming back the data:
    Y <- 2^Y+1
    assays(se)[[1]] <- Y
    return(se)
}




#' @title Within-sample scaling normalization for dual-guides screens.
#' 
#' @description Within-sample scaling normalization for dual-guides screens.
#' 
#' @param se A SummarizedExperiment object.
#' @param fun String specifying which function should be
#'     used for normalization. "median" by default.
#' @param type.field String specifying column name in 
#'     \code{ann} containing gRNA pairs class.
#'     "class" by default.
#' @param type.levels.dko Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the double
#'     knockout constructs. 
#' @param type.levels.sko1 Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the single
#'     knockout constructs at position 1.
#' @param type.levels.sko2 Character vector specifying which
#'     classes of gRNA pairs in \code{ann[[type.field]]}
#'     should be used for normalization of the single
#'     knockout constructs at position 2.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param gene.pair.field String specifying colum name in \code{rowData(se)}
#'     containing gene pair name.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @details The normalization makes the distributions of sKO1, sKO2 and dKO
#'     constructs comparable within a sample. The \code{type.levels.sko1},
#'     \code{type.levels.sko2} and \code{type.levels.dko} specify the pairs
#'     to be used for each class respectively to estimate the normalization
#'     scaling factors. This should be performed before running
#'     \code{normalizeDual}.
#' 
#' @return A SummarizedExperiment with normalized counts.
#' 
#' @examples
#' normalizeDualWithinSample(seDualExample)
#' 
#' @importFrom matrixStats colMedians
#' @export
normalizeDualWithinSample <- function(se,
                                      fun=c("median", "mean"),
                                      type.field="class",
                                      type.levels.dko=NULL,
                                      type.levels.sko1=NULL,
                                      type.levels.sko2=NULL,
                                      gene1.field="gene_symbol_1",
                                      gene2.field="gene_symbol_2",
                                      gene.pair.field="group",
                                      sko.control="neg"
){
    fun <- match.arg(fun)
    if (!type.field %in% colnames(rowData(se))){
        stop("type.field not found in rowData(se)")
    }
   

    # Extracting data:
    assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    Y <- assays(se)[[1]]
    Y <- log2(Y+1)

    sko1.indices <- .getSkoIndices(se=se,
                                   sko.position="first",
                                   sko.control=sko.control,
                                   gene1.field=gene1.field,
                                   gene2.field=gene2.field,
                                   gene.pair.field=gene.pair.field)
    sko2.indices <- .getSkoIndices(se=se,
                                   sko.position="second",
                                   sko.control=sko.control,
                                   gene1.field=gene1.field,
                                   gene2.field=gene2.field,
                                   gene.pair.field=gene.pair.field)
    dko.indices <- .getDkoIndices(se=se,
                                  sko.control=sko.control,
                                  gene1.field=gene1.field,
                                  gene2.field=gene2.field,
                                  gene.pair.field=gene.pair.field)
    wh.dko  <- dko.indices
    wh.sko1 <- sko1.indices
    wh.sko2 <- sko2.indices
    if (!is.null(type.levels.dko)){
        wh.dko <- which(rowData(se)[[type.field]]==type.levels.dko)
    }
    if (!is.null(type.levels.sko1)){
        wh.sko1 <- which(rowData(se)[[type.field]]==type.levels.sko1)
    }
    if (!is.null(type.levels.sko2)){
        wh.sko2 <- which(rowData(se)[[type.field]]==type.levels.sko2)
    }

    if (length(wh.dko)==0){
        stop("None of the type.levels.dko are found in the specified type.field column.")
    }
    if (length(wh.sko1)==0){
        stop("None of the type.levels.sko1 are found in the specified type.field column.")
    }
    if (length(wh.sko2)==0){
        stop("None of the type.levels.sko2 are found in the specified type.field column.")
    }
    if (fun=="median"){
        factors.dko  <- colMedians(as.matrix(Y)[wh.dko,,drop=FALSE], na.rm=TRUE)
        factors.sko1 <- colMedians(as.matrix(Y)[wh.sko1,,drop=FALSE], na.rm=TRUE)
        factors.sko2 <- colMedians(as.matrix(Y)[wh.sko2,,drop=FALSE], na.rm=TRUE)
    } else if (fun=="mean"){
        factors.dko  <- colMeans(as.matrix(Y)[wh.dko,,drop=FALSE], na.rm=TRUE)
        factors.sko1 <- colMeans(as.matrix(Y)[wh.sko1,,drop=FALSE], na.rm=TRUE)
        factors.sko2 <- colMeans(as.matrix(Y)[wh.sko2,,drop=FALSE], na.rm=TRUE)
    }
    factors <- rbind(factors.sko1, factors.sko2, factors.dko)
    factors <- apply(factors,2, function(x){
        x-median(x)
    })
    factors.sko1 <- factors[1,]
    factors.sko2 <- factors[2,]
    factors.dko  <- factors[3,]
    
    Y[sko1.indices,] <- sweep(Y[sko1.indices,,drop=FALSE],2,factors.sko1, "-")
    Y[sko2.indices,] <- sweep(Y[sko2.indices,,drop=FALSE],2,factors.sko2, "-")
    Y[dko.indices,]  <- sweep(Y[dko.indices,,drop=FALSE],2,factors.dko, "-")


    # Transforming back the data:
    Y <- 2^Y+1
    assays(se)[[1]] <- Y
    return(se)
}






# removeDkoEffects <- function(se,
#                              reference.level="Reference",
#                              condition.field="Group",
#                              replicate.field="Replicate",
#                              fun=c("median", "mean"),
#                              sko.control="neg",
#                              gene.pair.field="group",
#                              gene1.field="gene_symbol_1",
#                              gene2.field="gene_symbol_2"
# ){
#     fun <- match.arg(fun)
#     ratios <- .getLogRatiosForNormalization(se,
#                                             reference.level=reference.level,
#                                             reference.field=condition.field,
#                                             replicate.field=replicate.field)
#     lfc1 <- getSkoData(ratios,
#                        sko.position="first",
#                        return.matrix=TRUE,
#                        aggregate=TRUE)
#     lfc2 <- getSkoData(ratios,
#                        sko.position="second",
#                        return.matrix=TRUE,
#                        aggregate=TRUE)
#     lfc <- getDkoData(ratios,
#                       return.matrix=TRUE,
#                       aggregate=TRUE)
#     .getFirstGene <- function(names){
#         unlist(lapply(strsplit(names, split="_"), function(x) x[[1]]))
#     }

#     .getSecondGene <- function(names){
#         unlist(lapply(strsplit(names, split="_"), function(x) x[[2]]))
#     }

#     offsets <- vapply(seq_len(ncol(se)), function(i){
#         df <- data.frame(lfc=lfc[,i])
#         df$lfc1 <- lfc1[,i][.getFirstGene(rownames(lfc))]
#         df$lfc2 <- lfc2[,i][.getSecondGene(rownames(lfc))] 
#         if (fun=="median"){
#             offset <- median(df$lfc - (df$lfc1+df$lfc2), na.rm=TRUE)
#         } else if (fun=="mean"){
#             offset <- mean(df$lfc - (df$lfc1+df$lfc2), na.rm=TRUE)
#         }   
#     }, FUN.VALUE=0)
    
    
    
   
#     dko.indices <- .getDkoIndices(se,
#                                   sko.control=sko.control,
#                                   gene.pair.field=gene.pair.field,
#                                   gene1.field=gene1.field,
#                                   gene2.field=gene2.field)
#     #col.indices <- which(colData(se)[[condition.field]]!=reference.level)
#     Y <- assays(ratios)[[1]]
#     Y[dko.indices, ] <- sweep(Y[dko.indices,], 2, offsets, "-")
#     assays(ratios)[[1]] <- Y

#     se <- .getInverseLogRatiosForNormalization(se=se,
#                                                se.ratio=ratios,
#                                                reference.field=condition.field,
#                                                reference.level=reference.level,
#                                                replicate.field=replicate.field)
#     return(se)
# }



