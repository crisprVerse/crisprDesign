#' @title Add isoform-specific annotation to a \linkS4class{GuideSet} object
#' 
#' @description Add isoform-specific annotation to a
#'      \linkS4class{GuideSet} object.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param tx_id String specifiying Ensembl ID for the 
#'     isoform transcript of interested. E.g. "ENST00000311936".
#' @param ... Additional arguments, currently ignored.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @return A A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.with the 
#'     following added columns: \code{percentCDS},
#'     \code{percentCodingIsoforms}, and 
#'     \code{isCommonCodingExon}. The column values are
#'     specific to the transcript specified by \code{tx_id}.
#'     The \code{percentCDS} column indicates at what percentage
#'     of the coding sequence the gRNA is cutting. The
#'     column \code{percentCodingIsoforms} indicates the 
#'     percentage of coding isoforms that are targeted
#'     by the gRNA. The column \code{isCommonCodingExon}
#'     indicates whether or not the exon targetd by the
#'     gRNA is common to all isoforms for the gene.
#' 
#' @examples
#' data(guideSetExampleFullAnnotation)
#' tx_id <- "ENST00000538872"
#' gs <- guideSetExampleFullAnnotation
#' gs <- addIsoformAnnotation(gs, tx_id)
#' 
#' @rdname addIsoformAnnotation
#' @export
setMethod("addIsoformAnnotation",
          "GuideSet", 
          function(object,
                   tx_id
){
    if (!"geneAnnotation" %in% colnames(mcols(object))){
        stop("geneAnnotation not found in GuideSet. ",
             "See ?addGeneAnnotation to add a gene annotation first.")
    }
    cols <- c("percentCDS",
              "percentCodingIsoforms",
              "isCommonCodingExon",
              "pfam")
    geneAnn <- geneAnnotation(object, unlist=TRUE)

    # Let's extract gene id:
    wh <- match(tx_id, geneAnn$tx_id)
    if (length(wh)==0){
        stop("tx_id not found")
    }
    gene_id <- geneAnn[wh, "gene_id"]
    geneAnn <- geneAnn[geneAnn$gene_id==gene_id,,drop=FALSE]
    cols <- intersect(cols, colnames(geneAnn))

    # Extracting metadata columns:
    df <- S4Vectors::mcols(object)
    # Initializing
    for (col in cols){
        df[[col]] <- NA
    } 
    wh <- match(rownames(df), rownames(geneAnn))
    if ("percentCodingIsoforms" %in% cols){
        df$percentCodingIsoforms <- geneAnn[wh,"percentCodingIsoforms"]
    }
    if ("isCommonCodingExon" %in% cols){
        df$isCommonCodingExon <- geneAnn[wh,"isCommonCodingExon"]   
    }
    
    if (!is.na(tx_id)){
        cols <- c("percentCDS", "pfam")
        ann <- geneAnnotation(object,
                              tx_id=tx_id,
                              unlist=TRUE)
        cols <- intersect(cols, colnames(ann))
        if (nrow(ann)>0){
            ann <- ann[, cols, drop=FALSE]
            
            # Adding to object:
            wh <- match(rownames(df), rownames(ann))
            for (col in cols){
                df[[col]] <- ann[wh,col]
            }
        }
    }
    S4Vectors::mcols(object) <- df
    return(object)
})


#' @rdname addIsoformAnnotation
#' @export
setMethod("addIsoformAnnotation",
          "PairedGuideSet", 
          function(object,
                   tx_id
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addIsoformAnnotation(unifiedGuideSet,
                                            tx_id=tx_id)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    
    return(out)
})



#' @rdname addIsoformAnnotation
#' @export
setMethod("addIsoformAnnotation", "NULL", function(object){
    return(NULL)
})



