#' @title Return crisprBliss scores for dual-guide screens
#' 
#' @description Return crisprBliss scores for dual-guide screens.
#' 
#' @param se SummarizedExperiment object containing normalized counts.
#' @param group.field String specifying colum name in \code{colData(se)}
#'     containing sample group for comparisons.
#' @param group.reference String specifying the reference group 
#'    in \code{colData(se)[group.field]}.
#' @param group.later String specifying sample condition
#'    in \code{colData(se)[group.field]} to be contrasted to 
#'    \code{group.reference}.
#' @param replicate.field String specifying colum name in \code{colData(se)}
#'     containing sample replicate information.
#' @param gene.pairs data.frame with two columns specifying gene pairs
#'    for which to calculate crisprBliss scores. First column specifies
#'    the first gene, and the second column specified the second gene.
#'    See \code{getGenePairs} to obtain such a data.frame.
#' @param id1.field String specifying colum name in \code{rowData(se)} 
#'     containing ID for gRNA in position 1.
#' @param id2.field String specifying colum name in \code{rowData(se)} 
#'     containing ID for gRNA in position 2.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param robust Should robust regression be used? FALSE by default.
#' @param robust.method String specifing robust regression method
#'    if \code{robust=TRUE}. "Huber" by default.
#' @param control.type.field String specifying column name in 
#'     \code{ann} containing gRNA pairs class.
#'     "class" by default.
#' @param control.type.level String specifying which gRNA pairs
#'     in \code{ann[[control.type.field]]} should be used
#'     for estimating maximal dropout. 
#' @param quantile Numeric value specifying quantile to be used
#'     to calulcate maximal dropout. 0.05 by default.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param assay Numeric value specifying the index of the assay
#'     in \code{assays(se)} to be used. 1 by default.
#' @param verbose Should messaged be printed? TRUE by default.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @return A data.frame. 
#' 
#' @examples
#' out <- getCrisprBlissScores(seDualExample,
#'                             group.field="Group",
#'                             group.reference="Ref",
#'                             group.later="Day16")
#' 
#' @importFrom multcomp glht Chisqtest
#' @importFrom MASS rlm
#' @export
getCrisprBlissScores <- function(se,
                                 group.field="Group",
                                 group.reference="Reference",
                                 group.later=NULL,
                                 replicate.field="Replicate",
                                 gene.pairs=NULL,
                                 id1.field="ID_1",
                                 id2.field="ID_2",
                                 gene1.field="gene_symbol_1",
                                 gene2.field="gene_symbol_2",
                                 robust=FALSE,
                                 robust.method=c("huber","hampel", "bisquare"),
                                 control.type.field="class",
                                 control.type.level="eg_eg",
                                 quantile=0.01,
                                 sko.control="neg",
                                 assay=1,
                                 verbose=TRUE
){
    # Calculating offsets for dKO toxicity adjustment:
    wh <- which(colData(se)[[group.field]]==group.later)
    offsets <- .getDkoCuttingOffsets(se, 
                                     reference.level=group.reference,
                                     condition.field=group.field,
                                     replicate.field=replicate.field,
                                     fun="median")
    offsets <- offsets[wh]
    
    

    if (is.null(gene.pairs)){
        gene.pairs <- getGenePairs(se,
                                   gene1.field=gene1.field,
                                   gene2.field=gene2.field)
    }
    wh_top    <- which(colData(se)[[group.field]]==group.later)
    wh_bottom <- which(colData(se)[[group.field]]==group.reference)
    Y_top     <- log2(assays(se)[[assay]]+1)[,wh_top]
    Y_bottom  <- log2(assays(se)[[assay]]+1)[,wh_bottom]

    if (length(group.later)>1){
        stop("crisprBliss only works with group.later of length 1.")
    }
    ann <- as.data.frame(rowData(se))
    robust.method <- match.arg(robust.method)
    robust.method <- paste0("psi.", robust.method)

      #Estimate bliss.ceiling
    wh <- which(ann[[control.type.field]]==control.type.level)
    if (length(wh)==0){
        warning("Specified control.type.level is not found. Using",
                " target_target instead. You might want to change",
                " the default control.quantile.")
        wh <- which(ann[[control.type.field]]=="target_target")
        if (length(wh)==0){
            stop("cannot find target_target either.")
        }
    }

  
    #Estimate maximal dropout:
    dropout <- matrix(rowMeans(Y_top)-rowMeans(Y_bottom),ncol=1)
    maxDropout <- .getMaximalDropout(dropout, 
                                     ann=ann,
                                     control.type.field=control.type.field,
                                     control.type.level=control.type.level,
                                     quantile=quantile)
    if (verbose){
        message("Maximal dropout is estimated to be ", maxDropout)
    }

    


    .getPairStats <- function(ann,
                              gene1,
                              gene2,
                              robust=FALSE,
                              correct.dko.bias=FALSE,
                              offsets=NULL){

        .getDesignMatrix <- function(){
            oneSampleDesign <- .createOneSampleDesign()
            topDesign    <- .getMultiSamplesDesign(oneSampleDesign, ncol(Y_top))
            bottomDesign <- .getMultiSamplesDesign(oneSampleDesign, ncol(Y_bottom))
            fullDesign   <- .getFinalDesign(topDesign, bottomDesign)
            return(fullDesign)
        }

        # Getting single-knockout main terms:
        .createOneSampleDesign <- function(){
            guide1  <- as.character(info[[id1.field]])
            guide2  <- as.character(info[[id2.field]])
            design1 <- model.matrix(~guide1-1)
            design2 <- model.matrix(~guide2-1)
            design1 <- design1[, grepl(gene1, colnames(design1)),drop=FALSE]
            design2 <- design2[, grepl(gene2, colnames(design2)),drop=FALSE]
            design  <- cbind(design1, design2)
            colnames(design) <- gsub("guide1", "", colnames(design))
            colnames(design) <- gsub("guide2", "", colnames(design))

             # Getting double-knockout effect:
            x <- ifelse(c12,1,0)
            design <- cbind(design,x)
            colnames(design)[ncol(design)] <- "score"
            return(design)
        }
        

        .getMultiSamplesDesign <- function(design, n_samples){
            fullDesign <- design
            for (sample in 2:n_samples){
                fullDesign <- rbind(fullDesign, design)
            }
            fullDesign <- cbind(fullDesign, fullDesign)
            return(fullDesign)
        }

        .getFinalDesign <- function(topDesign,
                                    bottomDesign){
            n_coeffs <- ncol(topDesign)/2
            bottom_coeffs <- seq_len(n_coeffs)
            top_coeffs <- bottom_coeffs+n_coeffs
            bottomDesign[,top_coeffs] <- 0
            topDesign[,bottom_coeffs] <- 0
            fullDesign <- rbind(bottomDesign, topDesign)
            colnames(fullDesign)[bottom_coeffs] <- paste0("ref_",
                                                       colnames(fullDesign)[bottom_coeffs])
            colnames(fullDesign)[top_coeffs] <- paste0("tx_",
                                                       colnames(fullDesign)[top_coeffs])
            return(fullDesign)
        }


        .createLincomp <- function(fullDesign, bliss.factor){
            K <- matrix(0, nrow=5, ncol=ncol(fullDesign))
            colnames(K) <- colnames(fullDesign)
            rownames(K) <- c("gene1",
                             "gene2",
                             "synergy",
                             "dko",
                             "bliss")
            # sKO1:
            wh1_top    <- grepl(gene1, colnames(K)) & grepl("tx", colnames(K))
            wh1_bottom <- grepl(gene1, colnames(K)) & grepl("ref", colnames(K))
            K[1,wh1_top] <- 1/sum(wh1_top)
            K[1,wh1_bottom] <- -1/sum(wh1_bottom)
            # sKO2:
            wh2_top    <- grepl(gene2, colnames(K)) & grepl("tx", colnames(K))
            wh2_bottom <- grepl(gene2, colnames(K)) & grepl("ref", colnames(K))
            K[2,wh2_top] <- 1/sum(wh2_top)
            K[2,wh2_bottom] <- -1/sum(wh2_bottom)
            # synergy:
            wh3_top    <- grepl("score", colnames(K)) & grepl("tx", colnames(K))
            wh3_bottom <- grepl("score", colnames(K)) & grepl("ref", colnames(K))
            K[3,wh3_top]    <- 1/sum(wh3_top)
            K[3,wh3_bottom] <- -1/sum(wh3_bottom)

            # dKO:
            K[4,wh1_top] <- 1/sum(wh1_top)
            K[4,wh2_top] <- 1/sum(wh2_top)
            K[4,wh3_top] <- 1/sum(wh3_top)
            K[4,wh1_bottom] <- -1/sum(wh1_bottom)
            K[4,wh2_bottom] <- -1/sum(wh2_bottom)
            K[4,wh3_bottom] <- -1/sum(wh3_bottom)

            # Now creating Bliss combination:
            K[5,wh2_top] <- (1/sum(wh2_top))*(1-bliss.factor)
            K[5,wh3_top] <- 1/sum(wh3_top)
            K[5,wh2_bottom] <- -(1/sum(wh2_bottom))*(1-bliss.factor)
            K[5,wh3_bottom] <- -1/sum(wh3_bottom)
            return(K)
        }

       .getLfcVector <- function(){
            lfc1 <- as.vector(Y_bottom[rownames(info),])
            lfc2 <- as.vector(Y_top[rownames(info),])
            lfc <- c(lfc1,lfc2)
            return(lfc)
        }

        .getAdjustedLfcVector <- function(offsets, c12){
            lfc1 <- as.vector(Y_bottom[rownames(info),])
            lfc2 <- Y_top[rownames(info),,drop=FALSE]
            lfc2[c12,] <- sweep(lfc2[c12,],2,offsets,"-")
            lfc2 <- as.vector(lfc2)
            lfc <- c(lfc1,lfc2)
            return(lfc)
        }



       

        cond1 <- ann[[gene1.field]] %in% c(gene1, sko.control)
        cond2 <- ann[[gene2.field]] %in% c(gene2, sko.control)
        cond3 <- ann[[gene1.field]]==sko.control & ann[[gene2.field]]==sko.control
        info  <- ann[cond1 & cond2 &!cond3,]
        info  <- as.data.frame(info)
        patt <- paste0("_", sko.control)
        if (nrow(info)<=3 | sum(grepl(patt, rownames(info)))==1){
            return(rep(NA, 10))
        }
        patt1 <- paste0("_", sko.control, "$")
        patt2 <- paste0("^", sko.control, "_")
        c1  <- grepl(patt1, info$class)
        c2  <- grepl(patt2, info$class)
        c12 <- !c1 & !c2

        lfc1 <- Y_bottom[rownames(info),]
        lfc2 <- Y_top[rownames(info),]
        avg1  <- mean(lfc2[c1,,drop=FALSE], na.rm=TRUE)-mean(lfc1[c1,,drop=FALSE], na.rm=TRUE)
        avg2  <- mean(lfc2[c2,,drop=FALSE], na.rm=TRUE)-mean(lfc1[c2,,drop=FALSE], na.rm=TRUE)
        avg12 <- mean(lfc2[c12,,drop=FALSE], na.rm=TRUE)-mean(lfc1[c12,,drop=FALSE], na.rm=TRUE)
        bliss.factor <- .getBlissFactor(avg1, maxDropout)

        fullDesign <- .getDesignMatrix()
        if (correct.dko.bias){
            lfc <- .getAdjustedLfcVector(offsets=offsets,
                                         c12=c12)
        } else {
            lfc <- .getLfcVector()
        }
        

        if (!robust){
            model <- lm(lfc~fullDesign-1)
        } else {
            model <- MASS::rlm(lfc~fullDesign-1,
                               psi=robust.method)
        }
        
        # Create lincomp matrix:
        K <- .createLincomp(fullDesign, bliss.factor)
        coeffs <- rep(0, nrow(K))
        pvals  <- rep(0, nrow(K))
        for (k in seq_len(nrow(K))){
            temp <- summary(multcomp::glht(model, K[k,,drop=FALSE]),
                            test=Chisqtest())  
            coeffs[k] <- temp$test$coefficients
            pvals[k]  <- temp$test$pvalue
        }
        names(coeffs) <- names(pvals) <- rownames(K)
  
        out=c(score=coeffs["synergy"],
              pval=pvals["synergy"],
              bliss.score=coeffs["bliss"], 
              bliss.pval=pvals["bliss"],
              avg1, avg2, avg12,
              pval.sko1=pvals["gene1"],
              pval.sko2=pvals["gene2"],
              pval.dko=pvals["dko"])
        return(out)
    }

    if (verbose){
        cat("[getCrisprBlissScores] Calculating synergy scores. \n")
    }
    results <- lapply(seq_len(nrow(gene.pairs)), function(i){
        .getPairStats(ann,
                      gene.pairs[i,1],
                      gene.pairs[i,2],
                      robust=robust)
    })
    results <- do.call(rbind, results)
    rownames(results) <- rownames(gene.pairs)
    
    colnames(results) <- c("score.synergy.naive",
                           "pval.synergy.naive",
                           "score.synergy.bliss",
                           "pval.synergy.bliss",
                           "lfc.sko1",
                           "lfc.sko2",
                           "lfc.dko",
                           "pval.sko1",
                           "pval.sko2",
                           "pval.dko")
    results <- as.data.frame(results)
    pval.cols <- colnames(results)[grepl("pval", colnames(results))]
    for (col in pval.cols){
        fdr.col <- gsub("pval", "fdr", col)
        results[[fdr.col]] <- p.adjust(results[[col]], "fdr")    
    }


    if (verbose){
        cat("[getCrisprBlissScores] Calculating adjusted synergy scores. \n")
    }
    # Getting results adjusted for dKO toxicity:
    results.adjusted <- lapply(seq_len(nrow(gene.pairs)), function(i){
        .getPairStats(ann,
                      gene.pairs[i,1],
                      gene.pairs[i,2],
                      robust=robust,
                      correct.dko.bias=TRUE,
                      offsets=offsets)
    })
    results.adjusted <- do.call(rbind, results.adjusted)
    rownames(results.adjusted) <- rownames(gene.pairs)
    results.adjusted <- results.adjusted[,1:4]
    colnames(results.adjusted) <- c("score.synergy.naive.adjusted",
                                    "pval.synergy.naive.adjusted",
                                    "score.synergy.bliss.adjusted",
                                    "pval.synergy.bliss.adjusted")
    results.adjusted <- as.data.frame(results.adjusted)
    pval.cols <- colnames(results.adjusted)[grepl("pval", colnames(results.adjusted))]
    for (col in pval.cols){
        fdr.col <- gsub("pval", "fdr", col)
        results.adjusted[[fdr.col]] <- p.adjust(results.adjusted[[col]], "fdr")    
    }

    

    # Adding more scores:
    s  <- results[["score.synergy.naive"]]
    y1 <- results[["lfc.sko1"]]
    y2 <- results[["lfc.sko2"]]
    bliss <- results[["score.synergy.bliss"]]
    results$score.strong    <- abs(s)- pmax(abs(y1), abs(y2), na.rm=TRUE)
    results$score.lethality <- pmin(y1,y2) - (y1+y2+s)
    results$score.recovery  <- (y1+y2+s) - pmin(y1,y2)
    results$synergy.bliss.expected.dko <- s - bliss + y1 + y2
    results$targetPair <- rownames(results) 
    results$target1 <- .getFirstTarget(rownames(results))
    results$target2 <- .getSecondTarget(rownames(results))

    results <- cbind(results, results.adjusted)

    cols <- c("targetPair",
              "target1",
              "target2",
              "score.synergy.naive",
              "pval.synergy.naive",
              "fdr.synergy.naive",
              "score.synergy.bliss",
              "pval.synergy.bliss",
              "fdr.synergy.bliss",
              "synergy.bliss.expected.dko",
              "score.synergy.naive.adjusted",
              "pval.synergy.naive.adjusted",
              "score.synergy.bliss.adjusted",
              "pval.synergy.bliss.adjusted",
              "lfc.sko1",
              "pval.sko1",
              "fdr.sko1",
              "lfc.sko2",
              "pval.sko2",
              "fdr.sko2",
              "lfc.dko",
              "pval.dko",
              "fdr.dko",
              "score.strong",
              "score.lethality",
              "score.recovery")
    results <- results[,cols,drop=FALSE]
    return(results)
}




#' @importFrom SummarizedExperiment assays
#' @export
.getMaximalDropout <- function(ratios,
                               ann,
                               control.type.field="class",
                               control.type.level="eg_eg",
                               quantile=0.05
){
    if (is(ratios, "SummarizedExperiment")){
        ratios <- assays(ratios)[[1]]
    }
    if (!control.type.field %in% colnames(ann)){
        stop("control.type.field not in colnames(ann)")
    }
    wh_controls <- ann[[control.type.field]]==control.type.level
    if (length(wh_controls)==0){
        stop("control.type.level not found in ann[[control.type.field]]")
    }
    x <- rowMeans(ratios[wh_controls, , drop=FALSE], na.rm=TRUE)
    dropout <- quantile(x, quantile, na.rm=TRUE)
    return(dropout)
}



.getBlissFactor <- function(sko1,
                            maximalDropout,
                            force.finite=TRUE
){
    if (sko1>0){
        factor <- 1
    } else {
        factor <- 1-min(abs(sko1/maximalDropout),1)
    }
    if (!is.finite(factor) & force.finite){
        factor <- 1
    }
    return(factor)
}



.getDkoCuttingOffsets <- function(se,
                                  reference.level="Reference",
                                  condition.field="Condition",
                                  replicate.field="Replicate",
                                  fun=c("median", "mean")
){
    fun <- match.arg(fun)
    ratios <- .getLogRatiosForNormalization(se,
                                            reference.level=reference.level,
                                            reference.field=condition.field,
                                            replicate.field=replicate.field)
    lfc1 <- getSkoData(ratios,
                       sko.position="first",
                       return.matrix=TRUE,
                       aggregate=TRUE)
    lfc2 <- getSkoData(ratios,
                       sko.position="second",
                       return.matrix=TRUE,
                       aggregate=TRUE)
    lfc <- getDkoData(ratios,
                      return.matrix=TRUE,
                      aggregate=TRUE)
    .getFirstGene <- function(names){
        unlist(lapply(strsplit(names, split="_"), function(x) x[[1]]))
    }

    .getSecondGene <- function(names){
        unlist(lapply(strsplit(names, split="_"), function(x) x[[2]]))
    }


    offsets <- vapply(seq_len(ncol(se)), function(i){
        df <- data.frame(lfc=lfc[,i])
        df$lfc1 <- lfc1[,i][.getFirstGene(rownames(lfc))]
        df$lfc2 <- lfc2[,i][.getSecondGene(rownames(lfc))] 
        if (fun=="median"){
            offset <- median(df$lfc - (df$lfc1+df$lfc2), na.rm=TRUE)
        } else if (fun=="mean"){
            offset <- mean(df$lfc - (df$lfc1+df$lfc2), na.rm=TRUE)
        }   
    }, FUN.VALUE=0)
    names(offsets) <- colnames(se)
    return(offsets)
}









# #' @importFrom multcomp glht Chisqtest
# #' @importFrom MASS rlm
# #' @export
# runCrisprBliss <- function(Y,
#                            ann,
#                            gene.pairs,
#                            robust=FALSE,
#                            robust.method=c("huber","hampel", "bisquare"),
#                            control.type.field="class",
#                            control.type.level="eg_eg",
#                            control.quantile=0.05,
#                            sko.control="neg",
#                            verbose=TRUE
# ){
#     robust.method <- match.arg(robust.method)
#     robust.method <- paste0("psi.", robust.method)
  
#     #Estimate bliss.ceiling
#     wh <- which(ann[[control.type.field]]==control.type.level)
#     if (length(wh)==0){
#         warning("Specified control.type.level is not found. Using",
#                 " target_target instead. You might want to change",
#                 " the default control.quantile.")
#         wh <- which(ann[[control.type.field]]=="target_target")
#         if (length(wh)==0){
#             stop("cannot find target_target either.")
#         }
#     }
#     x <- rowMeans(Y[wh,,drop=FALSE], na.rm=TRUE)
#     bliss.ceiling <- quantile(x, control.quantile)
#     if (verbose){
#         message("Maximal dropout is estimated to be ", bliss.ceiling)
#     }

#     .getPairStats <- function(Y, ann, gene1, gene2, robust=FALSE){

#         .createLincomp <- function(fullDesign, bliss.factor){
#             K <- matrix(0, nrow=5, ncol=ncol(fullDesign))
#             colnames(K) <- colnames(fullDesign)
#             rownames(K) <- c("gene1", "gene2", "synergy", "dko", "bliss")
#             wh1 <- grepl(gene1, colnames(K))
#             wh2 <- grepl(gene2, colnames(K))
#             wh3 <- grepl("score", colnames(K))
#             K[1,wh1] <- 1/sum(wh1)
#             K[2,wh2] <- 1/sum(wh2)
#             K[3,wh3] <- 1/sum(wh3)
#             K[4,wh1] <- 1/sum(wh1)
#             K[4,wh2] <- 1/sum(wh2)
#             K[4,wh3] <- 1/sum(wh3)
#             # Now creating Bliss combination:
#             if (!is.finite(bliss.factor)){
#                 bliss.factor <- 1
#             }
#             K[5,wh2] <- (1/sum(wh2))*(1-bliss.factor)
#             K[5,wh3] <- 1/sum(wh3)
#             return(K)
#         }


#         cond1 <- ann$gene_symbol_1 %in% c(gene1, sko.control)
#         cond2 <- ann$gene_symbol_2 %in% c(gene2, sko.control)
#         cond3 <- ann$gene_symbol_1==sko.control & ann$gene_symbol_2==sko.control
#         info  <- ann[cond1 & cond2 &!cond3,]
#         info  <- as.data.frame(info)
#         patt <- paste0("_", sko.control)
#         if (nrow(info)<=3 | sum(grepl(patt, rownames(info)))==1){
#             return(rep(NA, 10))
#         }
#         patt1 <- paste0("_",sko.control, "$")
#         patt2 <- paste0("^",sko.control, "_")
#         c1  <- grepl(patt1, info$class)
#         c2  <- grepl(patt2, info$class)
#         c12 <- !c1 & !c2

#         lfc <- as.vector(Y[rownames(info),])
#         avg1   <- mean(lfc[c1], na.rm=TRUE)
#         avg2   <- mean(lfc[c2], na.rm=TRUE)
#         avg12  <- mean(lfc[c12], na.rm=TRUE)
#         bliss.factor <- 1-min(abs(avg1/bliss.ceiling),1)

#         # Getting single-knockout main terms:
#         guide1  <- as.character(info$ID_1)
#         guide2  <- as.character(info$ID_2)
#         design1 <- model.matrix(~guide1-1)
#         design2 <- model.matrix(~guide2-1)
#         design1 <- design1[, grepl(gene1, colnames(design1)),drop=FALSE]
#         design2 <- design2[, grepl(gene2, colnames(design2)),drop=FALSE]
#         design  <- cbind(design1, design2)
#         colnames(design) <- gsub("guide1", "", colnames(design))
#         colnames(design) <- gsub("guide2", "", colnames(design))

#         # Getting double-knockout effect:
#         x <- ifelse(c12,1,0)
#         design <- cbind(design,x)
#         colnames(design)[ncol(design)] <- "score"
#         nsamples <- ncol(Y)
       
#         fullDesign <- design
#         for (k in 2:nsamples){
#             fullDesign <- rbind(fullDesign,design)
#         }

       
#         # Fitting model and extracting effects:
#         #model  <- summary(lm(lfc~fullDesign-1))
#         if (!robust){
#             model <- lm(lfc~fullDesign-1)
#         } else {
#             model <- MASS::rlm(lfc~fullDesign-1, psi=robust.method)
#         }
        
#         # Create lincomp matrix:
#         K <- .createLincomp(fullDesign, bliss.factor)
#         coeffs <- rep(0, nrow(K))
#         pvals  <- rep(0, nrow(K))
#         for (k in seq_len(nrow(K))){
#             temp <- summary(multcomp::glht(model, K[k,,drop=FALSE]),
#                             test=Chisqtest())  
#             coeffs[k] <- temp$test$coefficients
#             pvals[k]  <- temp$test$pvalue
#         }
#         names(coeffs) <- names(pvals) <- rownames(K)
  
#         out=c(score=coeffs["synergy"],
#               pval=pvals["synergy"],
#               bliss.score=coeffs["bliss"], 
#               bliss.pval=pvals["bliss"],
#               avg1, avg2, avg12,
#               pval.sko1=pvals["gene1"],
#               pval.sko2=pvals["gene2"],
#               pval.dko=pvals["dko"])
#         return(out)
#     }

#     results <- lapply(seq_len(nrow(gene.pairs)), function(i){
#         #print(i)
#         .getPairStats(Y, ann, gene.pairs[i,1], gene.pairs[i,2], robust=robust)
#     })
#     results <- do.call(rbind, results)
#     rownames(results) <- rownames(gene.pairs)
#     colnames(results) <- c("synergy.lfc",
#                            "synergy.pval",
#                            "bliss.lfc",
#                            "bliss.pval",
#                            "sko1.lfc",
#                            "sko2.lfc",
#                            "dko.lfc",
#                            "sko1.pval",
#                            "sko2.pval",
#                            "dko.pval")
#     results <- as.data.frame(results)
#     results$synergy.fdr <- p.adjust(results$synergy.pval, "fdr")
#     results$bliss.fdr   <- p.adjust(results$bliss.pval, "fdr")
#     results$sko1.fdr    <- p.adjust(results$sko1.pval, "fdr")
#     results$sko2.fdr    <- p.adjust(results$sko2.pval, "fdr")
#     results$dko.fdr     <- p.adjust(results$dko.pval, "fdr")

#     # Adding more scores:
#     s  <- results$synergy.lfc
#     y1 <- results$sko1.lfc
#     y2 <- results$sko2.lfc
#     bliss <- results$bliss.lfc
#     results$score.strong    <- abs(s)- pmax(abs(y1), abs(y2), na.rm=TRUE)
#     results$score.lethality <- pmin(y1,y2) - (y1+y2+s)
#     results$score.recovery  <- (y1+y2+s) - pmin(y1,y2)
#     results$bliss.expected.dko <- s - bliss + y1 + y2
      

#     cols <- c("synergy.lfc", "synergy.pval", "synergy.fdr",
#              "bliss.lfc", "bliss.pval",
#              "bliss.fdr","bliss.expected.dko",
#              "sko1.lfc", "sko1.pval", "sko1.fdr",
#              "sko2.lfc", "sko2.pval", "sko2.fdr",
#              "dko.lfc",  "dko.pval", "dko.fdr",
#              "score.strong", "score.lethality", "score.recovery")
#     results <- results[,cols]
#     return(results)
# }

