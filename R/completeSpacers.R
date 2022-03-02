#' @title Get complete spacer information
#' @description These functions serve to "fill-in-the-blank" for
#'     spacers lacking information.
#' 
#' @param start Coordinate of the first nucleotide of the spacer sequences.
#'     Must be always less than \code{end}.
#' @param end Coordinate of the last nucleotide of the spacer sequence.
#'     Must be always greater than \code{start}.
#' @param chr The chromosome in which the protospacer sequence is located.
#' @param pam_site Coordinate of the first nucleotide of the PAM sequence.
#' @param strand Either "+" or "-".
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param bsgenome A \linkS4class{BSgenome} object.
#' @param spacerLen Spacer sequence length.
#'     If NULL, the information is obtained from \code{crisprNuclease}.
#' 
#' @return A numeric or character vector, depending on the function.
#' 
#' @return \code{getPAMSiteFromStartAndEnd} returns a numeric vector of the
#'     coordinate of the first nucleotide of the PAM sequences.
#' @return \code{getCutSiteFromPamSite} returns a numeric vector of coordinate
#'     of the cut site.
#' @return \code{getPAMSequence} returns a character vector of PAM sequences.
#' @return \code{getSpacerSequence} returns a character vector of
#'     spacer sequences.
#' 
#' @details Functions that return coordinates (\code{getPAMSite},
#'     \code{getCutSite}, \code{getSpacerRanges}) do not check whether
#'     coordinates exceed chromosomal lengths.
#' @details The start and end coordinates of a genomic range is
#'     strand-independent, and always obeys \code{start <= end}.
#'
#' @examples
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' dat <- data.frame(chr='chr4', start=1642343, strand='+')
#' dat$pam_site <- getPAMSiteFromStartAndEnd(start=dat$start,
#'                                           strand=dat$strand)
#' dat$cut_site <- getCutSiteFromPamSite(pam_site=dat$pam_site,
#'                                       strand=dat$strand)
#' dat$pam <- getPAMSequence(chr=dat$chr,
#'                           pam_site=dat$pam_site,
#'                           strand=dat$strand,
#'                           bsgenome=bsgenome)
#' dat$spacer <- getSpacerSequence(chr=dat$chr,
#'                                 pam_site=dat$pam_site,
#'                                 strand=dat$strand,
#'                                 bsgenome=bsgenome)
#' }
#' @name completeSpacers
NULL









#' @rdname completeSpacers
#' @export
#' @importFrom crisprBase pamLength
#' @importFrom crisprBase pamSide
#' @importFrom crisprBase spacerLength
getPAMSiteFromStartAndEnd <- function(start=NULL,
                                      end=NULL,
                                      strand,
                                      crisprNuclease=NULL,
                                      spacerLen=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pamlen   <- pamLength(crisprNuclease)
    pamside  <- pamSide(crisprNuclease)
    # check input
    .checkStartEnd(start,end)
    strand <- .validateStrand(strand)

    in_lengths <- c(length(start), length(end), length(strand))
    in_lengths <- unique(in_lengths[in_lengths > 0])
    if (length(in_lengths) > 1){
        stop('start and end must be the same length as strand or null.')
    }
    if (is.null(spacerLen)){
        spacerLen <- spacerLength(crisprNuclease)
    }
    # initialize PAM site coordinates
    pam_site <- rep(0, in_lengths)
    pos <- which(strand=='+')        # positive/forward strand
    neg <- which(strand=='-')        # negative/reverse strand
    # get start/end if either is null
    if (is.null(start)){
        start <- end - spacerLen + 1
    }
    if (is.null(end)){
        end <- start + spacerLen - 1
    }
    # get coordinates by strand
    if (pamside=="3prime"){
        pam_site[pos] <- end[pos] + 1
        pam_site[neg] <- start[neg] - 1
    } else {
        pam_site[pos] <- start[pos] - pamlen
        pam_site[neg] <- end[neg] + pamlen
    }
    return(pam_site)
}






#' @rdname completeSpacers
#' @export
getCutSiteFromPamSite <- function(pam_site,
                                  strand,
                                  crisprNuclease=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    # pam_site <- .validatePamSite(pam_site)
    strand <- .validateStrand(strand)
    stopifnot(length(pam_site)==length(strand))
    cut_offset <- .getDefaultCutOffset(crisprNuclease)
    cut_site <- pam_site
    cut_site[strand=='+'] <- pam_site[strand=='+'] + cut_offset
    cut_site[strand=='-'] <- pam_site[strand=='-'] - cut_offset
    return(cut_site)
}





#' @rdname completeSpacers
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome getSeq
#' @export
getPAMSequence <- function(chr,
                           pam_site,
                           strand, 
                           crisprNuclease=NULL,
                           bsgenome=NULL
){
    # handle inputs
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pamlen   <- pamLength(crisprNuclease)
    pamside  <- pamSide(crisprNuclease)
    chr      <- unlist(chr)
    pam_site <- unlist(pam_site)
    strand   <- unlist(strand)
    if (length(unique(vapply(c(chr,
                               pam_site,
                               strand), length, FUN.VALUE=1)))!=1){
        stop('chr, pam_site, and strand must be the same length.')
    }
    pam_site <- .validatePamSite(pam_site)
    strand   <- .validateStrand(strand)


    bad_chr <- !chr %in% names(seqlengths(bsgenome))
    if (sum(bad_chr) > 0){
        stop('chr name not recognized: ', paste(chr[bad_chr], collapse=', '))
    }
    chr_lens <- seqlengths(bsgenome)[chr]
    if (!is.numeric(pam_site) || sum(pam_site < 1) > 0 || 
        sum(pam_site > chr_lens) > 0){
        stop('invalid PAM site.')
    }

    # get range of PAM
    start <- end <- pam_site #Initializing
    start[strand=='+'] <- pam_site[strand=='+'] 
    end[strand=='+']   <- pam_site[strand=='+'] + pamlen - 1
    start[strand=='-'] <- pam_site[strand=='-'] - pamlen + 1
    end[strand=='-']   <- pam_site[strand=='-'] 
   
    pam <- getSeq(bsgenome,
                  chr,
                  start=start,
                  end=end,
                  strand=strand)
    pam <- as.character(pam)
    return(pam)
}



#' @importFrom crisprBase pamLength
getPAMSequence_customSeq <- function(custom_seq,
                                     pam_site,
                                     strand,
                                     crisprNuclease=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    start <- end <- pam_site
    pamlen <- pamLength(crisprNuclease)
    start[strand == "-"] <- pam_site[strand == "-"] - pamlen + 1
    end[strand == "+"]   <- pam_site[strand == "+"] + pamlen - 1 
    pam <- vapply(seq_along(pam_site), function(i){
        seq <- substr(custom_seq, start[i], end[i])
        if (strand[[i]]=="-"){
            seq <- .revComp(seq)
        }
        return(seq)
    }, FUN.VALUE="a")
    pam <- as.character(pam)
    return(pam)
}





#' @rdname completeSpacers
#' @export
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom BSgenome getSeq
getSpacerSequence <- function(chr,
                              pam_site,
                              strand,
                              crisprNuclease=NULL,
                              bsgenome=NULL,
                              spacerLen=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pamlen   <- pamLength(crisprNuclease)
    pamside  <- pamSide(crisprNuclease)
    # handle input
    if (unique(vapply(c(chr,
                        pam_site,
                        strand), length, FUN.VALUE=1))!=1){
        stop('chr, pam_site, and strand must have the same length.')
    }
    pam_site <- .validatePamSite(pam_site)
    strand   <- .validateStrand(strand)


    bad_chr <- !chr %in% names(seqlengths(bsgenome))
    if (sum(bad_chr) > 0){
        stop('chr name not recognized: ', paste(chr[bad_chr], collapse=', '))
    }
    chr_lens <- seqlengths(bsgenome)[chr]
    if (sum(pam_site < 1) > 0 || sum(pam_site > chr_lens) > 0){
        stop('invalid PAM site.')
    }
    if (is.null(spacerLen)){
        spacerLen <- spacerLength(crisprNuclease)
    } 
    spacerLen <- .validateSpacerLength(spacerLen)
    # get spacer range limits
    if (pamside=="3prime"){
        start <- pam_site - spacerLen
        end   <- pam_site - 1
        start[strand=='-'] <- pam_site[strand=='-'] + 1
        end[strand=='-']   <- pam_site[strand=='-'] + spacerLen
    } else {
        start <- pam_site + pamlen
        end   <- pam_site + (pamlen + spacerLen) -1
        start[strand=='-'] <- pam_site[strand=='-'] - (pamlen + spacerLen) + 1
        end[strand=='-']   <- pam_site[strand=='-'] - pamlen
    }
    # get spacer sequence
    seqs <- getSeq(bsgenome, chr, start=start, end=end, strand=strand)
    seqs <- as.character(seqs)
    return(seqs)
}








#' @importFrom S4Vectors metadata DataFrame
#' @importFrom BiocGenerics start end start<- end<- strand
#' @importFrom crisprBase pamSide pamLength spacerLength
convertToProtospacerGRanges <- function(guideSet){
    guideSet <- .validateGuideSet(guideSet)
    crisprNuclease <- crisprNuclease(guideSet)
    pamSide  <- pamSide(crisprNuclease)
    spacer_len <- spacerLength(crisprNuclease)
    pam_len <- pamLength(crisprNuclease)

    # Get protospacer GRanges:
    r <-  as.character(BiocGenerics::strand(guideSet))=='-'
  
    if (pamSide=='3prime'){
        start    <- pamSites(guideSet) - spacer_len
        end      <- pamSites(guideSet) + pam_len -1
        start[r] <- pamSites(guideSet)[r] - pam_len +1
        end[r]   <- pamSites(guideSet)[r] + spacer_len
    } else {
        start    <- pamSites(guideSet)
        end      <- pamSites(guideSet) + spacer_len + pam_len -1 
        start[r] <- pamSites(guideSet)[r] - spacer_len - pam_len + 1
        end[r]   <- pamSites(guideSet)[r]
    }
    gr.new <- guideSet
    BiocGenerics::start(gr.new) <- start
    BiocGenerics::end(gr.new)   <- end
    return(gr.new)
}






#' @importFrom IRanges trim
#' @importFrom BSgenome getSeq
.getExtendedSequences <- function(guideSet,
                                  start,
                                  end
){
    guideSet <- .validateGuideSet(guideSet)
    
    gr <- guideSet
    wh_neg <- which(as.character(strand(gr))=="-")
    # The order of resizing IRanges matters
    # to presever the validity of a positive width.
    if (start>0 & end>0){
        end(gr)   <- end(guideSet)+end
        start(gr) <- start(guideSet)+start
        start(gr)[wh_neg] <- start(guideSet)[wh_neg]-end
        end(gr)[wh_neg]   <- end(guideSet)[wh_neg]-start
    } else {
        start(gr) <- start(guideSet)+start
        end(gr)   <- end(guideSet)+end
        end(gr)[wh_neg]   <- end(guideSet)[wh_neg]-start
        start(gr)[wh_neg] <- start(guideSet)[wh_neg]-end
    }

    gr <- GenomicRanges::trim(gr) #Taking care of invalid values
    if (targetOrigin(guideSet)=="customSequences"){
        seqs <- getSeq(customSequences(guideSet),gr)
    } else {
        seqs <- getSeq(bsgenome(guideSet), gr)
    }
    seqs <- as.character(seqs)

    #Making sure the sequences are not out of bound:
    len = end-start+1 # Expected length
    seqs[seqs==""] <- NA
    seqs[nchar(seqs)<len] <- NA
    return(seqs)
}






.validatePamSite <- function(pam_site){
    stopifnot(is.numeric(unlist(pam_site)))
    if (sum(pam_site < 1) > 0 || sum(pam_site %% 1 != 0) > 0){
        stop('pam_site must be a positive integer.')
    }
    return(pam_site)
}

.validateStrand <- function(strand){
    stopifnot(is.character(strand))
    if (sum(!strand %in% c('+', '-')) > 0){
        stop('strand must contain either "+" or "-" values.')
    }
    return(strand)
}

.validateSpacerLength <- function(len){
    stopifnot(is.numeric(len))
    if (length(len) > 1){
        stop('len accepts a single value only.')
    }
    if (len < 1 || len %% 1 != 0){
        stop('len must be a positive integer.')
    }
    return(len)
}



#' @importFrom crisprBase cutSites
.getDefaultCutOffset <- function(crisprNuclease){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    offset <- cutSites(crisprNuclease,
                       middle=TRUE)
    return(offset)
}


.checkStartEnd <- function(start, end){
    if (is.null(start) && is.null(end)){
        stop('start and/or end must be provided.')
    }
    if (!is.null(start)){
        stopifnot(is.numeric(start))
        if (sum(start <= 0) > 0 || sum(start %% 1 != 0) > 0){
            stop('start must contain positive integers only.')
        }
    }
    if (!is.null(end)){
        stopifnot(is.numeric(end))
        if (sum(end <= 0) > 0 || sum(end %% 1 != 0) > 0){
            stop('end must contain positive integers only.')
        }
    }
}





# #' @rdname completeSpacers
# #' @importFrom IRanges IRanges
# #' @export
# getSpacerRanges <- function(pam_site,
#                             strand,
#                             nuclease=NULL,
#                             len=NULL
# ){
#     nuclease <- .validateNuclease(nuclease)
#     pamlen   <- pamLength(nuclease)
#     pamside  <- pamSide(nuclease)
#     # check input
#     strand <- .validateStrand(strand)
#     pam_site <- .validatePamSite(pam_site)
#     stopifnot(length(pam_site)==length(strand))
#     if (is.null(len)){
#         len <- spacerLength(nuclease)
#     } 
#     len <- .validateSpacerLength(len)
#     # get start and end coordinates
#     start <- end <- pam_site
#     if (pamside == "3prime"){
#         start[strand=='+'] <- pam_site[strand=='+'] - len
#         start[strand=='-'] <- pam_site[strand=='-'] + 1
#         end[strand=='+']   <- pam_site[strand=='+'] - 1
#         end[strand=='-']   <- pam_site[strand=='-'] + len
#     } else {
#         start[strand=='+'] <- pam_site[strand=='+'] + pamlen
#         start[strand=='-'] <- pam_site[strand=='-'] - pamlen - len + 1
#         end[strand=='+']   <- pam_site[strand=='+'] + pamlen + len - 1
#         end[strand=='-']   <- pam_site[strand=='-'] - pamlen
#     }
#     ir <- IRanges(start=start, end=end)
#     return(ir)
# }





# #' @rdname completeSpacers
# #' @importFrom IRanges IRanges
# #' @export
# getProtospacerRanges <- function(pam_site,
#                                  strand,
#                                  nuclease=NULL,
#                                  len=NULL
# ){
#     nuclease <- .validateNuclease(nuclease)
#     pamlen   <- pamLength(nuclease)
#     pamside  <- pamSide(nuclease)
#     # check input
#     strand <- .validateStrand(strand)
#     pam_site <- .validatePamSite(pam_site)
#     stopifnot(length(pam_site)==length(strand))
#     if (is.null(len)){
#         len <- spacerLength(nuclease)
#     } 
#     len <- .validateSpacerLength(len)
#     # get start and end coordinates
#     start <- end <- pam_site
#     if (pamside == "3prime"){
#         start[strand=='+'] <- pam_site[strand=='+'] - len
#         start[strand=='-'] <- pam_site[strand=='-'] - pamlen + 1
#         end[strand=='+']   <- pam_site[strand=='+'] + pamlen - 1 
#         end[strand=='-']   <- pam_site[strand=='-'] + len
#     } else {
#         start[strand=='+'] <- pam_site[strand=='+'] 
#         start[strand=='-'] <- pam_site[strand=='-'] - (pamlen + len) + 1
#         end[strand=='+']   <- pam_site[strand=='+'] + (pamlen + len) - 1
#         end[strand=='-']   <- pam_site[strand=='-'] 
#     }
#     ir <- IRanges(start=start, end=end)
#     return(ir)
# }

