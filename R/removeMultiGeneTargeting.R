#' @title Remove gRNAs targeting secondary targets
#' 
#' @description Remove gRNAs targeting secondary targets
#' 
#' @param guideSet A \linkS4class{GuideSet} object. 
#' @param geneID String specifying gene ID of the main gene target.
#' @param geneColumn Column in \code{geneAnnotation(guideSet)} specifying
#'     gene IDs
#' @param ignoreGenesWithoutSymbols Should gene without gene symbols
#'     be ignored when removing co-targeting gRNAs?
#' @param ignoreReadthroughs Should readthrough genes be ignored when 
#'     removing co-targeting gRNAs?
#' 
#' @details The protospacer target sequence of gRNAs can be located in 
#'    overlapping genes, and this function allows users to filter out such 
#'    gRNAs. This ensures remaining gRNAs are targeting only one gene.
#' 
#' @return A \linkS4class{GuideSet} object with gRNAs targeting multiple
#'     targets removed. 
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @rdname removeSpacersWithSecondaryTargets
removeSpacersWithSecondaryTargets <- function(guideSet,
                                              geneID,
                                              geneColumn="gene_id",
                                              ignoreGenesWithoutSymbols=TRUE,
                                              ignoreReadthroughs=TRUE
){
    stopifnot(.hasGeneAnnotation(guideSet))
    geneAnn <- geneAnnotation(guideSet)
    if (!geneColumn %in% colnames(geneAnn)){
        stop("Column specified by geneColumn not found ",
             "in geneAnnotation(guideSet).")
    }
    geneAnn <- geneAnn[geneAnn[[geneColumn]]!=geneID,,drop=FALSE]
   

    if (ignoreGenesWithoutSymbols & nrow(geneAnn)!=0){
        if ("gene_symbol" %in% colnames(geneAnn)){
            geneAnn <- geneAnn[!is.na(geneAnn$gene_symbol),,drop=FALSE]
            geneAnn <- geneAnn[which(geneAnn$gene_symbol!=""),,drop=FALSE]
        }
    }

    if (ignoreReadthroughs & nrow(geneAnn)!=0){
        if ("gene_symbol" %in% colnames(geneAnn)){
            geneAnn <- geneAnn[!grepl("-", geneAnn$gene_symbol),,drop=FALSE]
        }
    }
    if (nrow(geneAnn)==0){
        out <- guideSet 
    } else {
        badGuides <- unique(rownames(geneAnn))
        out <- guideSet[!names(guideSet) %in% badGuides]
    }
    return(out)
}
