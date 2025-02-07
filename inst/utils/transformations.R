# condition.field="Condition"
# reference.level="Reference"
# replicate.field="Replicate"
# #replicate.field=NULL
# seLogRatios <- getLogRatios(seExample, 
#   condition.field="Condition",
#   reference.level="Reference",
#   replicate.field="Replicate"
# )
#' @title Get log2 fold changes of samples with respect to a reference group
#' @description Computes the log2 fold changes of guides for all samples with
#' respect to the reference timepoint or condition and returns a
#' SummarizedExperiment with those values.
#' @param se A SummarizedExperiment object.
#' @param condition.field String specifying the column of \code{colData(se)}
#'     containing the sample group information.
#' @param reference.level String specifying the reference level of the column
#'     named by \code{condition.field}.
#' @param replicate.field String specifying the column of \code{colData(se)}
#'     containing the replicate sample information.
#' @param use.normalized Should normalized counts be used? TRUE by default.
#' @return A SummarizedExperiment object with log2 fold change data with
#' respect to the reference. By definition, reference samples will have values
#' of \code{0} for all guides. The metadata for the SummarizedExperiment will
#' also have "log_ratio" in the \code{transform} list, which is used by other
#' functions to determine whether the input SummarizedExperiment has (or does
#' not have) the necessary data transformations.
#' @examples 
#' seLog <- getLogRatios(seExample, 
#'  condition.field="Condition",
#'  reference.level="Reference",
#'  replicate.field="Replicate",
#'  use.normalized=FALSE
#' )
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
getLogRatios <- function(se,
                         condition.field=NULL,
                         reference.level=NULL,
                         replicate.field=NULL,
                         use.normalized=TRUE
){

    .checkConditionField(se, condition.field)
    .checkReferenceLevel(se, condition.field, reference.level)
    .checkReplicateField(se, replicate.field)

    group <- colData(se)[,condition.field]
    if (is.null(replicate.field)){
        replicate <- rep("Rep1", ncol(se))   
    } else {
        replicate <- colData(se)[,replicate.field]
    }
    replicate_levels <- unique(replicate)

    
    n_samples <- ncol(se)
    refSamples <- lapply(seq_len(n_samples), function(sample){
        sampleSpecificReplicate <- replicate[sample]
        wh_replicate <- replicate==sampleSpecificReplicate
        wh_group <- group==reference.level
        wh <- which(wh_replicate & wh_group)
        if (length(wh)!=1){
            message("At least one sample does not have an obvious reference
                    sample, so ignoring replicate.")
            wh <- which(group==reference.level)
        }
        colnames(se)[wh]
    })
     
    
    if (use.normalized){
        if (!"logcounts" %in% names(assays(se))){
            stop("logcounts assays not found and is required for use.normalized=TRUE ")
        } else {
            logCounts <- as.matrix(assays(se)[["logcounts"]])
        }
    } else {
        logCounts <- log2(as.matrix(assays(se)[[1]])+1)
    }

    # get log ratio matrix
    logRatios <- matrix(
        rep(NA, prod(dim(logCounts))),
            nrow=nrow(logCounts),
            dimnames=list(rownames(logCounts), colnames(logCounts)))
        # take log ratios
    for (sample in seq_len(n_samples)){
        denominator <- rowMeans(logCounts[,refSamples[[sample]],drop=FALSE])
        logRatios[,sample] <- logCounts[,sample] - denominator
    }

    # add transformation label to metadata
    metadata(se)$transform <- c(metadata(se)$transform, 'logRatio')

    # return SE object
    seLog <- SummarizedExperiment(assays=SimpleList(counts=logRatios),
                                  rowData=rowData(se),
                                  colData=colData(se),
                                  metadata=metadata(se))
    return(seLog)
}




.checkConditionField <- function(se,
                                 condition.field=NULL
){
    if (is.null(condition.field)){
        stop("condition.field cannot be NULL")
    }
    if (!condition.field %in% colnames(colData(se))){
        stop("condition.field is not found in colnames(se")
    }
    invisible(TRUE)
}

.checkReferenceLevel <- function(se,
                                 condition.field=NULL,
                                 reference.level=NULL
){
    .checkConditionField(se, condition.field)
    if (is.null(reference.level)){
        stop("reference.level cannot be NULL")
    }
    if (!reference.level %in% colData(se)[,condition.field]){
        stop("reference.field is not found in colData(se)[,condition.field]")
    }
    invisible(TRUE)
}


.checkReplicateField <- function(se,
                                 replicate.field
){
    if (is.null(replicate.field)){
        return(invisible(TRUE))
    } else if (!replicate.field %in% colnames(colData(se))){
        stop("replicate.field is not found in colnames(se")
    }
    invisible(TRUE)
}





#' @title Aggregate counts or log-ratios at the gene-level
#' @description Computes a summary statistic of counts or log-ratios across
#' all guides for all genes and returns a SummarizedExperiment object with 
#' gene-level aggregated values.
#' @param se A SummarizedExperiment object containing guide-level data
#'     (counts or log-ratios).
#' @param guess.log.transform Should appropriate log-transformations prior
#'     to aggregation be guessed from the data? FALSE by default. 
#' @param aggregate.field Column name in rowData specifying the
#'     grouping levels.
#' @param aggregate.ntcs Should non-targeting controls (NTCs) be aggregated?
#'     FALSE dy default.
#' @param fun Either "median" (default) or "mean". The summary statistic used 
#' to calculate gene-level data.
#' @param verbose Should messaged be printed to console? TRUE by default.
#' @return A SummarizedExperiment object with gene-level data. 
#'     The metadata for the SummarizedExperiment will also include
#'     "aggregate" in the \code{transform} list. 
#' @export
#' @importFrom stringr str_locate
#' @importFrom dplyr bind_cols
#' @importFrom magrittr %>%
#' @import S4Vectors
aggregateGuides <- function(se,
                            guess.log.transform=FALSE,
                            aggregate.field="group",
                            fun=c('median', 'mean'),
                            aggregate.ntcs=FALSE,
                            verbose=TRUE
){
    if (!guess.log.transform & verbose){
        cat(paste0("[aggregateGuides] Log transformations will not be inferred. ",
            "Users have to double-check that data have been transformed appropriately ",
            "before aggregation. \n"))
    }
    fun <- match.arg(fun)  
    ann <- rowData(se)
    if (!aggregate.field %in% colnames(ann)){
        stop("gene.field must be a column name in rowData")
    }
    grouping.var <- ann[[aggregate.field]]
    #If we want to keep NTCs separate: 
    if (!aggregate.ntcs){
        wh <- which(grouping.var=="NTC")
        if (length(wh)>0){
            if ("ID" %in% colnames(ann)){
                grouping.var[wh] <- ann[wh, "ID"]
            } else if ("id" %in% colnames(ann)){
                grouping.var[wh] <- ann[wh, "id"]
            } 
        }
        ann[[aggregate.field]] <- grouping.var
    }

    # Aggregating data:
    Ys <- assays(se)
    if (guess.log.transform){
        isLog  <- sapply(Ys, function(x){
            max(x, na.rm=TRUE)<30
        })
    } else {
        isLog <- rep(FALSE, length(Ys))
    }
    Ys[isLog] <- lapply(Ys[isLog], function(Y){
        return(2^Y-1)
    })
    Ys <- lapply(Ys,
                 .aggregateByRow,
                 fact=grouping.var,
                 fun=fun)
    Ys[isLog] <- lapply(Ys[isLog], function(Y){
        log2(Y+1)
    })
    Ys <- lapply(Ys, function(Y){
        colnames(Y) <- colnames(se)
        return(Y)
    })
    names(Ys) <- names(assays(se))

    
    # Modify feature data:
    cols <- c(aggregate.field,"group",
              "gene.symbol",
              "gene.id",
              "gene_symbol",
              "gene_id")
    cols <- intersect(cols, colnames(ann))
    ann_new <- ann[, cols, drop = FALSE]
    ann_new <- ann_new[!duplicated(ann_new[[aggregate.field]]),,drop=FALSE]
    rownames(ann_new) <- ann_new[[aggregate.field]]
    ann_new <- ann_new[rownames(Ys[[1]]),,drop=FALSE]
  
    # add transformation label to metadata
    metadata(se)$transform <- c(metadata(se)$transform, 'aggregate')
    metadata(se)$aggregate_function <- c(metadata(se)$aggregate_function, fun)
    metadata(se)$aggregate_ntc <- aggregate.ntcs
  
    # build aggregate summarized experiment  
    se_gene <- SummarizedExperiment(assays=Ys,
                                    rowData=ann_new,
                                    colData=colData(se),
                                    metadata=metadata(se))
    return(se_gene)
}


