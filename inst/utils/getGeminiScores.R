#' @title Return Gemini scores for dual-guide screens
#' 
#' @description Return Gemini scores for dual-guide screens.
#' 
#' @param se SummarizedExperiment object containing normalized counts.
#' @param group.field String specifying colum name in \code{colData(se)}
#'     containing sample group for comparisons.
#' @param group.reference String specifying the reference group 
#'    in \code{colData(se)[group.field]}.
#' @param groups.later Character vector specifying sample conditions
#'    in \code{colData(se)[group.field]} to be contrasted to 
#'    \code{group.reference}.
#' @param replicate.field String specifying colum name in \code{colData(se)}
#'     containing sample replicate information.
#' @param id1.field String specifying colum name in \code{rowData(se)} 
#'     containing ID for gRNA in position 1.
#' @param id2.field String specifying colum name in \code{rowData(se)} 
#'     containing ID for gRNA in position 2.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param n_iterations Number of iterations to use for Gemini
#'     convergence. See Gemini package.
#' @param force_results Should convergence be forced when running
#'     Gemini? TRUE by default. See Gemini package for more information.
#' @param assay Numeric value specifying the index of the assay
#'     in \code{assays(se)} to be used. 1 by default.
#' @param verbose Should messaged be printed? TRUE by default.
#' @param n_cores Integer specifing number of cores to be used
#'     in the Gemini computations. 12 by default.
#' @param offset Numeric value indicating an offset in the internal
#'     calculations of log fold changes in Gemini. 1 by default.
#' @examples
#' \dontrun{
#' se <- normalizeDual(seDualExample,
#'                     type.field="class",
#'                     type.levels="or_or")
#' 
#' out <- getGeminiScores(se,
#'                        group.field="Group",
#'                        group.reference="Ref",
#'                        groups.later=c("Day16","Day20_DMSO"))
#' }
#' @author Jean-Philippe Fortin
#' 
#' @return A list of \code{data.frame} objects. Each list element
#'     corresponds to a comparison specified by \code{groups.later}
#'     and \code{group.reference}. 
#'     Each data.frame has 3 columns corresponding the 3 scores
#'     calculated by GEMINI: strong, lethality, and recovery scores.
#' 
#' @importFrom gemini gemini_calculate_lfc
#' @importFrom gemini gemini_initialize
#' @importFrom gemini gemini_inference
#' @importFrom gemini gemini_score
#' @importFrom stats lm model.matrix quantile
#' @export
getGeminiScores <- function(se,
                            group.field="Group",
                            group.reference="Reference",
                            groups.later=NULL,
                            replicate.field="Replicate",
                            id1.field="ID_1",
                            id2.field="ID_2",
                            gene1.field="gene_symbol_1",
                            gene2.field="gene_symbol_2",
                            sko.control="neg",
                            n_iterations=100,
                            force_results=TRUE,
                            assay=1,
                            verbose=FALSE,
                            n_cores=12,
                            offset=1
){
    pattern_split <- ";"
    if (!group.field %in% colnames(colData(se))){
        stop("group.field should be part of colData(se).")
    }
    groups <- union(group.reference, groups.later)
    se <- se[,colData(se)[[group.field]] %in% groups]
    
    # Creating inputs:
    input <- .gemini_prepareInputs(se,
                                   group.field=group.field,
                                   replicate.field=replicate.field,
                                   pattern.split=pattern_split,
                                   assay=assay)
    counts.matrix               <- input$counts.matrix
    sample.replicate.annotation <- input$sample.replicate.annotation
    guide.annotation            <- input$guide.annotation
    colnames(counts.matrix) <- sample.replicate.annotation$colname
    counts.matrix <- as.matrix(counts.matrix)
    gene.column.names  <- c("gene_1", "gene_2")
    sample.column.name <- "samplename"
    etps <- which(sample.replicate.annotation$samplename==group.reference)
    ETP.column <- sample.replicate.annotation$colname[etps]

    g.input <- gemini::gemini_create_input(counts.matrix=counts.matrix,
                                           sample.replicate.annotation=sample.replicate.annotation,
                                           guide.annotation=guide.annotation,
                                           ETP.column=ETP.column,
                                           LTP.column=NULL,
                                           sample.column.name="sample.column.name",
                                           gene.column.names=gene.column.names,
                                           verbose=verbose)

    g.input.norm <- gemini::gemini_calculate_lfc(g.input, 
                                                 normalize=FALSE, 
                                                 CONSTANT=offset)
    g.model.norm <- gemini::gemini_initialize(Input=g.input.norm, 
                                              nc_gene=sko.control, 
                                              pattern_join=pattern_split,
                                              pattern_split=pattern_split, 
                                              cores=n_cores,
                                              verbose=verbose)
    g.inference <- gemini::gemini_inference(g.model.norm, 
                                            cores=n_cores, 
                                            n_iterations=n_iterations,
                                            force_results=force_results,
                                            verbose=verbose)
    g.score  <- gemini::gemini_score(Model=g.inference)

    # Reformatting
    comps <- colnames(g.score[[1]])
    comps <- paste0(comps,"_vs_",group.reference)
    ncomps <- ncol(g.score[[1]])
    out <- vector(length=ncomps, mode="list")
    names(out) <- comps
    pairs <- rownames(g.score[[1]])
    pairs <- gsub(";", "_", pairs)
    pairs <- .fixGeminiPairNames(pairs,
                                 ann=rowData(se),
                                 gene1.field=gene1.field,
                                 gene2.field=gene2.field)
    nscores <- length(g.score)
    scores <- names(g.score)
    for (comp in seq_len(ncomps)){
        out[[comp]] <- data.frame(targetPair=pairs)
        out[[comp]][["target1"]] <- .getFirstTarget(pairs)
        out[[comp]][["target2"]] <- .getSecondTarget(pairs)
        rownames(out[[comp]]) <- pairs
        for (score in seq_len(nscores)){
            score_name <- paste0("score.", scores[score])
            out[[comp]][, score_name] <- g.score[[scores[score]]][,comp]
        }
        colnames(out[[comp]]) <- gsub("_", '.', colnames(out[[comp]]))
        colnames(out[[comp]]) <- gsub("sensitive.", "", colnames(out[[comp]]))
    }
    return(out)
}

.getFirstTarget <- function(pairs){
    pairs <- strsplit(pairs, split="_")
    out <- lapply(pairs, function(x) x[[1]])
    unlist(out)
}

.getSecondTarget <- function(pairs){
    pairs <- strsplit(pairs, split="_")
    out <- lapply(pairs, function(x) x[[2]])
    unlist(out)
}


# Because GEMINI sorts pair by gene names -_-
.fixGeminiPairNames <- function(pairs,
                                ann,
                                gene1.field,
                                gene2.field
){
    df <- ann[, c(gene1.field, gene2.field)]
    df <- df[!duplicated(df),]
    df$pair <- paste0(df$gene_symbol_1, "_", df$gene_symbol_2)
    df$pair_rev <- paste0(df$gene_symbol_2, "_", df$gene_symbol_1)
    newPairs <- pairs
    wh <- which(!pairs %in% df$pair)
    newPairs[wh] <- df$pair[match(pairs[wh], df$pair_rev)]
    return(newPairs)
}




#Generate Input for Gemini
#se summarized experiment of raw count data
#assay assay number in se object
#replicate.field replicate field in rowData(se)
#group.field group field in rowData(se)
.gemini_prepareInputs <- function(se,
                                  id1.field="ID_1",
                                  id2.field="ID_2",
                                  gene1.field="gene_symbol_1",
                                  gene2.field="gene_symbol_2",
                                  assay=1,
                                  pattern.split=";",
                                  replicate.field="Replicate",
                                  group.field="Group"
){
    counts <- .gemini_prepareCounts(se,
                                    id1.field=id1.field,
                                    id2.field=id2.field,
                                    assay=assay,
                                    pattern.split=pattern.split)
    guide.annotation=.gemini_prepareGuideAnnotation(se,
                                                    id1.field=id1.field,
                                                    id2.field=id2.field,
                                                    gene1.field=gene1.field,
                                                    gene2.field=gene2.field,
                                                    pattern.split=pattern.split)
    sample.replicate.annotation=.gemini_prepareReplicateAnnotation(se,
                                                                   replicate.field=replicate.field,
                                                                   group.field=group.field)
    colnames(counts) <- sample.replicate.annotation[,"colname"]
    out <- list(counts.matrix=counts,
                guide.annotation=guide.annotation,
                sample.replicate.annotation=sample.replicate.annotation)
    return(out)
}






# Generate counts object for gemini
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
.gemini_prepareCounts <- function(se,
                                  id1.field="ID_1",
                                  id2.field="ID_2",
                                  assay=1,
                                  pattern.split=";"){
    ann <- rowData(se)
    df  <- assays(se)[[assay]]
    if (!id1.field %in% colnames(ann)){
        stop(paste0(id1.field, " not found in the rowData(se)"))
    }
    if (!id2.field %in% colnames(ann)){
        stop(paste0(id2.field, " not found in the rowData(se)"))
    }
    rownames(df) <- paste0(ann[[id1.field]], pattern.split, ann[[id2.field]])
    return(df)
}

# Generate guide annotation object for gemini
#' @importFrom SummarizedExperiment rowData
.gemini_prepareGuideAnnotation <- function(se,
                                           id1.field="ID_1",
                                           id2.field="ID_2",
                                           gene1.field="gene_symbol_1",
                                           gene2.field="gene_symbol_2",
                                           pattern.split=";"
){
    ann <- rowData(se)
    if (!id1.field %in% colnames(ann)){
        stop(paste0(id1.field, " not found in the rowData(se)"))
    }
    if (!id2.field %in% colnames(ann)){
        stop(paste0(id2.field, " not found in the rowData(se)"))
    }
    if (!gene1.field %in% colnames(ann)){
        stop(paste0(gene1.field, " not found in the rowData(se)"))
    }
    if (!gene2.field %in% colnames(ann)){
        stop(paste0(gene2.field, " not found in the rowData(se)"))
    }
    df <- data.frame(rowname=paste0(ann[[id1.field]], pattern.split, ann[[id2.field]]),
                     gene_1=as.character(ann[[gene1.field]]),
                     gene_2=as.character(ann[[gene2.field]]),
                     gene_1.guide=as.character(ann[[id1.field]]),
                     gene_2.guide=as.character(ann[[id2.field]]),
                     stringsAsFactors=FALSE)
    return(df)
}

# Generate sample replicate info object for gemini
#' @importFrom SummarizedExperiment colData
.gemini_prepareReplicateAnnotation <- function(se,
                                               replicate.field="Replicate",
                                               group.field="Group"){
    if (!replicate.field %in% colnames(colData(se))){
        stop(paste0(replicate.field, ", is not found in colData(se)"))
    }
    if (!group.field %in% colnames(colData(se))){
        stop(paste0(group.field, ", is not found in colData(se)"))
    }
    samplename <- colData(se)[ ,group.field]
    replicate  <- colData(se)[ ,replicate.field]
    colname    <- paste0(samplename, "_",replicate)
    df <- data.frame(colname=colname,
                     samplename=samplename,
                     replicate=replicate,
                     stringsAsFactors=FALSE)
    return(df)
}


# Get gene pairs for gemini's gemini_score function.
# Class should be NEG_NEG for the score function 
# (negative control pairs used for FDR calculation, etc.)
# pattern.split should be the same as used in gemini_initialize
# and in the annotation object created for gemini
#' @importFrom SummarizedExperiment rowData
gemini_prepareNCPairs <- function(se,
                                  type.field="class",
                                  type.level="or_or",
                                  gene1.field="gene_symbol_1",
                                  gene2.field="gene_symbol_2",
                                  pattern.split=";"
){
    ann <- rowData(se) 
    if (!type.field %in% colnames(ann)){
        stop(paste0(type.field, " is not a column found in rowData(se)"))
    }
    if (!gene1.field %in% colnames(ann)){
        stop(paste0(gene1.field, " not found in the rowData(se)"))
    }
    if (!gene2.field %in% colnames(ann)){
        stop(paste0(gene2.field, " not found in the rowData(se)"))
    }
    df    <- ann[ann[[type.field]]==type.level,,drop=FALSE]
    pairs <- paste0(df[[gene1.field]],pattern.split, df[[gene2.field]])
    pairs <- as.character(unique(pairs))
    return(pairs)
}




.mergeGeminiResults <- function(results){
    
    ID <- rownames(results[[1]])
    results <- lapply(results, function(x){
        x$ID <- NULL
        x
    })

    for (k in 1:length(results)){
        comp <- names(results)[[k]]
        colnames(results[[k]]) <- paste0(colnames(results[[k]]), 
                                         "_",comp)
    }
    names(results) <- NULL
    results <- do.call(cbind, results,)
    colnames(results) <- gsub("score", "gemini_", colnames(results))
    return(results)
}


