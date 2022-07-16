#' @author Jean-Philippe Fortin
#' 
#' @export 
flattenGuideSet <- function(guideSet){
    primaryTable <- .getPrimaryTable(guideSet)
    cols <- c("alignments",
              "geneAnnotation",
              "tssAnnotation",
              "enzymeAnnotation", "snps")
    cols <- intersect(cols, colnames(mcols(guideSet)))
    secondaryTables <- lapply(cols, function(col){
        .getSecondaryTable(guideSet, col)
    })
    names(secondaryTables) <- cols
    tables <- c(list(primary=primaryTable),
                secondaryTables)
    return(tables)
}



.getPrimaryTable <- function(guideSet){
    tab1 <- .getIrangesTable(guideSet)
    tab2 <- .getMcolsTable_flat(guideSet)
    tab <- cbind(tab1, tab2)
    tab <- .safeFormatColumns(tab)
    if (!"ID" %in% colnames(tab)){
        tab$ID <- rownames(tab)
    }
    rownames(tab) <- NULL
    tab <- .putColumnFirst("chr", tab)
    tab <- .putColumnFirst("ID", tab)
    return(tab)
}

#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
.getIrangesTable <- function(guideSet){
    out <- data.frame(chr=as.character(GenomeInfoDb::seqnames(guideSet)))
    out$start <- as.integer(BiocGenerics::start(guideSet))
    out$end   <- as.integer(BiocGenerics::end(guideSet))
    out$strand <- as.character(BiocGenerics::strand(guideSet))
    return(out)
}




#' @importFrom S4Vectors mcols
.getMcolsTable_flat <- function(guideSet){
    meta <- S4Vectors::mcols(guideSet)
    meta <- as.data.frame(meta)
    coltypes <-.getDFColtypes(meta)
    wh <- which(coltypes=="DNAStringSet")
    for (k in seq_along(wh)){
        meta[,wh[k]] <- as.character(meta[,wh[k]])
    }
    coltypes <-.getDFColtypes(meta)
    wh <- which(coltypes %in% .coltypes_flat)
    if (length(wh)>0){
        meta <- meta[,wh]    
    } else {
        meta <- NULL
    }
    return(meta)
}



#' @importFrom methods is
.getSecondaryTable <- function(guideSet,
                               colname=NULL
){

    if (colname=="alignments"){
        out <- crisprDesign::alignments(guideSet,
                                        unlist=TRUE)
    } else if (colname=="geneAnnotation"){
        out <- crisprDesign::geneAnnotation(guideSet,
                                            unlist=TRUE)
    } else if (colname=="tssAnnotation"){
        out <- crisprDesign::tssAnnotation(guideSet,
                                           unlist=TRUE)
    } else  if (colname=="snps"){
        out <- crisprDesign::snps(guideSet,
                                  unlist=TRUE)
    } else if (colname=="enzymeAnnotation"){
        out <- crisprDesign::enzymeAnnotation(guideSet,
                                              unlist=TRUE)
    } else {
        stop("colname not found in colnames(mcols(guideset)).")
    }
    if (is(out, "GRanges")){
        out$ID <- names(out)
        out <- .getPrimaryTable(out)
    } else {
        out$ID <- rownames(out)
        rownames(out) <- NULL
        out <- as.data.frame(out)
        if ("seqnames" %in% colnames(out)){
            out$chr <- out$seqnames
            out$seqnames <- NULL
        }
    }
    out <- .safeFormatColumns(out)
    if ("chr" %in% colnames(out)){
        out <- .putColumnFirst("chr", out)    
    }
    if ("ID" %in% colnames(out)){
        out <- .putColumnFirst("ID", out)    
    }
    return(out)
}


.putColumnFirst <- function(col, df){
    cols <- c(col, setdiff(colnames(df), col))
    df <- df[,cols, drop=FALSE]
    return(df)
}



.safeFormatColumns <- function(df){
    cols <- intersect(.coltypes_integer, colnames(df))
    if (length(cols)>0){
        for (k in seq_along(cols)){
            df[,cols[k]] <- as.integer(df[,cols[k]])
        }
    }
    return(df)
}



.getDFColtypes <- function(df){
    types <- vapply(seq_len(ncol(df)), function(i){
        class(df[,i])
    }, FUN.VALUE="a")
    return(types)
}


cols_aln <- paste0("n", 0:20)
cols_aln_c <- paste0(cols_aln, "_c")
cols_aln_p <- paste0(cols_aln, "_p")
cols_aln <- c(cols_aln, cols_aln_c, cols_aln_p)
.coltypes_integer <- c("start", "end", "pam_site", "cut_site",
                       "n_mismatches", "intergenic_distance",
                       "length","aminoAcidIndex",
                       "downstreamATG", "nIsoforms",
                       "totalIsoforms","nCodingIsoforms",
                       "tss_pos", "anchor_site", "dist_to_tss")
.coltypes_integer <- c(.coltypes_integer, cols_aln)
.coltypes_flat <- c("numeric", "logical", "integer", "character")

