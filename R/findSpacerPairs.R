#' @title Find pairs of CRISPR gRNA spacers from a pair of genomic regions.
#' 
#' @description Returns all possible, valid gRNA sequences for a given CRISPR
#'    nuclease from either a \linkS4class{GRanges} object or a set of
#'    sequence(s) contained in either a \linkS4class{DNAStringSet},
#'    \linkS4class{DNAString} or character vector of genomic sequences.
#' 
#' @param x1 Either a \linkS4class{GRanges}, a \linkS4class{DNAStringSet}, or a
#'     \linkS4class{DNAString} object, or a character vector of genomic
#'     sequences. This specifies the sequence space from which gRNAs in
#'     position 1 of the pairs will be designed. 
#' @param x2 Either a \linkS4class{GRanges}, a \linkS4class{DNAStringSet}, or a
#'     \linkS4class{DNAString} object, or a character vector of genomic
#'     sequences. This specifies the sequence space from which gRNAs in
#'     position 2 of the pairs will be designed. 
#' @param sortWithinPair Should gRNAs be sorted by chr and position
#'     within a pair? TRUE by default.
#' @param pamOrientation String specifying a constraint on the PAM orientation  
#'     of the pairs. Should be either "all" (default), "out" (for the so-called
#'     PAM-out orientation) or "in" (for PAM-in orientation). 
#' @param bsgenome A \linkS4class{BSgenome} object from which to extract
#'     sequences if \code{x} is a \linkS4class{GRanges} object.
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param canonical Whether to return only guide sequences having canonical
#'     PAM sequences. If TRUE (default), only PAM sequences with the highest
#'     weights stored in the \code{crisprNuclease} object will be considered.
#' @param spacer_len Length of spacers to return, if different from the
#'     default length specified by \code{crisprNuclease}.
#' @param both_strands Whether to consider both strands in search for 
#'     protospacer sequences. \code{TRUE} by default.
#' @param strict_overlap Whether to only include gRNAs that cut in the input
#'     range, as given by \code{cut_site} (\code{TRUE}) or to include all
#'     gRNAs that share any overlap with the input range (\code{FALSE}).
#'     \code{TRUE} by default. Ignored when \code{x} is not a
#'     \linkS4class{GRanges} object.
#' @param remove_ambiguities Whether to remove spacer sequences that contain
#'     ambiguous nucleotides (not explicily \code{A}, \code{C}, \code{G}, or
#'     \code{T}). TRUE by default.
#' 
#' @return A \linkS4class{PairedGuideSet} object.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @details This function returns a \linkS4class{PairedGuideSet} object 
#'     that stores gRNA pairs targeting the two genomic regions provided as 
#'     input. The gRNAs in position 1 target the first genomic region,
#'     and the gRNAs in position 2 target the second genomic region.
#'
#'     This function can be used for the following scenarios:
#'      
#'     1. Designing pairs of gRNAs targeting different genes, for instance
#'        for dual-promoter Cas9 systems, or polycystronic Cas12a constructs.
#'        This can also be used to target a given gene with multiple gRNAs
#'        for improved efficacy (for instance CRISPRa and CRISPRi)
#' 
#'     2. Designing pairs of gRNAs for double nicking systems such as Cas9 D10A.
#'     
#'     See vignette for more examples. 
#'     
#' @examples 
#' 
#' library(GenomicRanges)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' library(crisprBase)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' 
#' # Region 1:
#' gr1 <- GRanges(c("chr12"),
#'                IRanges(start=22224014, end=22225007))
#' 
#' # Region 2:
#' gr2 <- GRanges(c("chr12"),
#'                IRanges(start=23224014, end=23225007))
#' 
#' # Pairs targeting the same region:
#' pairs <- findSpacerPairs(gr1, gr1, bsgenome=bsgenome)
#' 
#' # Pairs targeting two regions:
#' # The gRNA in position targets gr1
#' # and the gRNA in position 2 targets gr2
#' pairs <- findSpacerPairs(gr1, gr2, bsgenome=bsgenome)
#' 
#' @seealso \code{\link{findSpacers}} to find unpaired spacer sequences, and the 
#'     \code{PairedGuideSet} object documentation to understand the
#'     output of \code{findSpacerPairs}.
#' 
#' @importFrom BiocGenerics sort
#' @export
findSpacerPairs <- function(x1,
                            x2,
                            sortWithinPair=TRUE,
                            pamOrientation=c("all", "out", "in"),
                            crisprNuclease=NULL,
                            bsgenome=NULL,
                            canonical=TRUE,
                            both_strands=TRUE,
                            spacer_len=NULL,
                            strict_overlap=TRUE,
                            remove_ambiguities=TRUE
){
    pamOrientation <- match.arg(pamOrientation)
    gs1 <- findSpacers(x1, 
                       crisprNuclease=crisprNuclease,
                       bsgenome=bsgenome,
                       canonical=canonical,
                       both_strands=both_strands,
                       spacer_len=spacer_len,
                       strict_overlap=strict_overlap,
                       remove_ambiguities=remove_ambiguities)
    gs2 <- findSpacers(x2, 
                       crisprNuclease=crisprNuclease,
                       bsgenome=bsgenome,
                       canonical=canonical,
                       both_strands=both_strands,
                       spacer_len=spacer_len,
                       strict_overlap=strict_overlap,
                       remove_ambiguities=remove_ambiguities)
    indices <- .getValidPairIndices(gs1,
                                    gs2,
                                    sortWithinPair=sortWithinPair)
        
    gs1 <- gs1[indices[,1]]
    gs2 <- gs2[indices[,2]]
    pgs <- PairedGuideSet(gs1, gs2)
    if (pamOrientation=="in"){
        pgs <- pgs[pamOrientation(pgs)=="in"]
    } else if (pamOrientation=="out"){
        pgs <- pgs[pamOrientation(pgs)=="out"]
    }
    return(pgs)
}


.orderingByCutSite <- function(gs){
    gs[order(cutSites(gs))]
}



.getValidPairIndices <- function(gs1,
                                 gs2,
                                 sortWithinPair
){
    id1 <- seq_along(gs1)
    id2 <- seq_along(gs2)
    df <- expand.grid(index1=id1, index2=id2)
    df$chr1 <- as.character(seqnames(gs1))[df$index1]
    df$chr2 <- as.character(seqnames(gs2))[df$index2]
    cuts1 <- cutSites(gs1)
    cuts2 <- cutSites(gs2)
    df$cut1 <- cuts1[df$index1]
    df$cut2 <- cuts2[df$index2]

    if (sortWithinPair){

        # Sorting by coordinates
        sameChr <- which(df$chr1==df$chr2)
        bad <- which(df[sameChr,]$cut1>=df[sameChr,]$cut2)
        if (length(bad)>0){
            df <- df[-bad,,drop=FALSE]
        }

        # Sorting by chromosome
        interChr <- which(df$chr1!=df$chr2)
        levels <- unique(c(df$chr1, df$chr2))
        df$chr1 <- factor(df$chr1, levels=levels)
        df$chr2 <- factor(df$chr2, levels=levels)
        df$chr1 <- as.numeric(df$chr1)
        df$chr2 <- as.numeric(df$chr2)
        bad <- which(df[interChr,]$chr1>df[interChr,]$chr2)
        if (length(bad)>0){
            df <- df[-bad,,drop=FALSE]
        }

    }
    return(df)
}






