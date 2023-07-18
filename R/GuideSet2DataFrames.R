#' @title Create a list of annotation tables from a GuideSet object
#' 
#' @description Create a list of annotation tables from a GuideSet object
#' 
#' @param guideSet A GuideSet object
#' @param useSpacerCoordinates Should the spacer coordinates be used
#'     as start and end coordinates? TRUE by default. If FALSE,
#'     the PAM site coordinate is used for both start and end. 
#' @param primaryOnly Should only the primary table (on-targets) be returned?
#'     FALSE by default.
#' 
#' @return A simple list of tables containing annotations derived from a 
#'     GuideSet object. The first table ("primary") is always available,
#'     while the other tables will be only available when the annotations
#'     were added to the GuideSet object. 
#' 
#' \itemize{
#' \item \code{primary} Primary table containing genomic coordinates and 
#'      sequence information of the gRNA sequences. Also contains on-target
#'      and off-target scores when available. 
#' \item \code{alignments} Table of on- and off-target alignments.
#' \item \code{geneAnnotation} Gene context annotation table.
#' \item \code{tssAnnotation} TSS context annotation table.
#' \item \code{enzymeAnnotation} Boolean table indicating whether or not 
#'     recognition motifs of restriction enzymes are found. 
#' \item \code{snps} SNP annotation table (human only).
#' }
#' 
#' @author Jean-Philippe Fortin
#' 
#' 
#' @export 
#' 
flattenGuideSet <- function(guideSet,
                            useSpacerCoordinates=TRUE,
                            primaryOnly=FALSE
){
    .Deprecated("GuideSet2DataFrames",
                msg="Deprecated soon. Use GuideSet2DataFrames instead.")
}




#' @title Create a list of annotation tables from a GuideSet object
#' 
#' @description Create a list of annotation tables from a GuideSet object
#' 
#' @param guideSet A GuideSet object
#' @param useSpacerCoordinates Should the spacer coordinates be used
#'     as start and end coordinates? TRUE by default. If FALSE,
#'     the PAM site coordinate is used for both start and end. 
#' @param primaryOnly Should only the primary table (on-targets) be returned?
#'     FALSE by default.
#' 
#' @return A simple list of tables containing annotations derived from a 
#'     GuideSet object. The first table ("primary") is always available,
#'     while the other tables will be only available when the annotations
#'     were added to the GuideSet object. 
#' 
#' \itemize{
#' \item \code{primary} Primary table containing genomic coordinates and 
#'      sequence information of the gRNA sequences. Also contains on-target
#'      and off-target scores when available. 
#' \item \code{alignments} Table of on- and off-target alignments.
#' \item \code{geneAnnotation} Gene context annotation table.
#' \item \code{tssAnnotation} TSS context annotation table.
#' \item \code{enzymeAnnotation} Boolean table indicating whether or not 
#'     recognition motifs of restriction enzymes are found. 
#' \item \code{snps} SNP annotation table (human only).
#' }
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' 
#' @export 
#' 
#' @examples
#' 
#' data(guideSetExampleFullAnnotation)
#' tables <- GuideSet2DataFrames(guideSetExampleFullAnnotation)
#' 
GuideSet2DataFrames <- function(guideSet,
                                useSpacerCoordinates=TRUE,
                                primaryOnly=FALSE
){
    primaryTable <- .getPrimaryTable(guideSet,
                                     useSpacerCoordinates=useSpacerCoordinates,
                                     addSpacer=TRUE)
    
    out <- list(primary=primaryTable)
    if (!primaryOnly){
        cols <- c("alignments",
                  "geneAnnotation",
                  "tssAnnotation",
                  "enzymeAnnotation",
                  "snps",
                  "txTable",
                  "exonTable",
                  "editedAlleles")
        cols <- intersect(cols, colnames(mcols(guideSet)))
        secondaryTables <- lapply(cols, function(col){
            .getSecondaryTable(guideSet,
                               colname=col,
                               useSpacerCoordinates=useSpacerCoordinates)
        })
        names(secondaryTables) <- cols
        out <- c(out, secondaryTables)
    }
    return(out)
}





.getPrimaryTable <- function(guideSet,
                             useSpacerCoordinates=TRUE,
                             nuclease=NULL,
                             addSpacer=TRUE
){
    if (is.null(nuclease)){
        nuclease <- crisprNuclease(guideSet)
    }
    tab1 <- .getIrangesTable(guideSet,
                             useSpacerCoordinates=useSpacerCoordinates,
                             nuclease=nuclease)
    tab2 <- .getMcolsTable_flat(guideSet)
    tab <- cbind(tab1, tab2)
    tab <- .safeFormatColumns(tab)
    if (!"ID" %in% colnames(tab)){
        tab[["ID"]] <- rownames(tab)
    }
    rownames(tab) <- NULL
    
    tab <- .putColumnFirst("chr", tab)
    if (is(guideSet, "GuideSet")){
        tab <- .putColumnFirst("protospacer", tab)        
        if (addSpacer){
            tab[["spacer"]] <- as.character(spacers(guideSet))
            tab <- .putColumnFirst("spacer", tab)
        }
    }
    tab <- .putColumnFirst("ID", tab)
    return(tab)
}


#' @importFrom GenomeInfoDb seqnames genome
#' @importFrom BiocGenerics start end strand
#' @importFrom crisprBase getProtospacerRanges
.getIrangesTable <- function(guideSet,
                             useSpacerCoordinates=TRUE,
                             nuclease=NULL
){
    out <- data.frame(chr=as.character(GenomeInfoDb::seqnames(guideSet)))
    if (useSpacerCoordinates){
        genome <- GenomeInfoDb::genome(guideSet)
        if (any(genome != "ntc")){
            validSeqnames <- names(genome[genome != "ntc"])
            protospacers <- as.vector(GenomeInfoDb::seqnames(guideSet)) %in%
                validSeqnames
            if (any(protospacers)){
                gr <- guideSet[protospacers]
                protospacerRanges <- crisprBase::getProtospacerRanges(
                    gr=gr,
                    nuclease=nuclease)
                guideSet[protospacers] <- protospacerRanges
            }
        }
    }
    out$start <- as.integer(BiocGenerics::start(guideSet))
    out$end   <- as.integer(BiocGenerics::end(guideSet))
    out$strand <- as.character(BiocGenerics::strand(guideSet))
    return(out)
}



#' @importFrom S4Vectors mcols
.getMcolsTable_flat <- function(guideSet){
    meta <- S4Vectors::mcols(guideSet)
    coltypes <- .getDFColtypes(meta)
    wh <- which(coltypes == "DNAStringSet")
    for (k in seq_along(wh)){
        meta[, wh[k]] <- as.character(meta[, wh[k]])
    }
    coltypes <- .getDFColtypes(meta)
    wh <- which(coltypes %in% .coltypes_flat)
    if (length(wh) > 0){
        meta <- meta[, wh, drop=FALSE]    
    } else {
        meta <- NULL
    }
    meta <- as.data.frame(meta)
    return(meta)
}



#' @importFrom methods is
.getSecondaryTable <- function(guideSet,
                               colname=NULL,
                               useSpacerCoordinates=TRUE
){
    nuclease <- crisprNuclease(guideSet)
    out <- switch(colname,
                  "alignments"=alignments(guideSet, unlist=TRUE),
                  "geneAnnotation"=geneAnnotation(guideSet, unlist=TRUE),
                  "tssAnnotation"=tssAnnotation(guideSet, unlist=TRUE),
                  "snps"=snps(guideSet, unlist=TRUE),
                  "enzymeAnnotation"=enzymeAnnotation(guideSet, unlist=TRUE),
                  "txTable"=txTable(guideSet, unlist=TRUE),
                  "exonTable"=exonTable(guideSet, unlist=TRUE),
                  "editedAlleles"=editedAlleles(guideSet, unlist=TRUE),
                  NULL)
    stopifnot("colname not found in colnames(mcols(guideset))." = {
        !is.null(out)
    })
    if (methods::is(out, "GRanges")){
        out$ID <- names(out)
        out <- .getPrimaryTable(out,
                                useSpacerCoordinates=useSpacerCoordinates,
                                nuclease=nuclease,
                                addSpacer=FALSE)
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
    df <- df[, cols, drop=FALSE]
    return(df)
}



.safeFormatColumns <- function(df){
    cols <- intersect(.coltypes_integer, colnames(df))
    if (length(cols) > 0){
        for (k in seq_along(cols)){
            df[, cols[k]] <- as.integer(df[, cols[k]])
        }
    }
    return(df)
}



.getDFColtypes <- function(df){
    types <- vapply(seq_len(ncol(df)), function(i){
        class(df[, i])
    }, FUN.VALUE=character(1))
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
