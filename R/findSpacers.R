#' @title Find CRISPR gRNA spacer sequences in a set of DNA sequences.
#' @description Returns all possible, valid gRNA sequences from either a
#'    \linkS4class{GRanges} object or a character vector of genomic sequences,
#'    for a given CRISPR nuclease.
#' 
#' @param x Either a \linkS4class{GRanges} object or a character vector of
#'     genomic sequences.
#' @param bsgenome \linkS4class{BSgenome} object from which to extract
#'     sequences if a \linkS4class{GRanges} object is provided as input. 
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param canonical Whether to return only guide sequences having canonical
#'     PAM sequences.
#' @param spacer_len Length of spacer to be used if different from
#'     default length specified in \code{crisprNuclease}.
#' @param both_strands Should both strands be considered to search for 
#'     protospacer sequences?
#'     TRUE by default.
#' @param cut_offset Distance in nucleotides between \code{pam_site}
#'     and \code{cut_site} if different from default offset specified by
#'     \code{crisprNuclease}.
#' @param strict_overlap If \code{TRUE} (default), only include guides
#'     having the cut site in the input range. Ignored when input is a
#'     character vector of genomic sequences.
#' @param collapse For input containing multiple regions, whether to return
#'     only unique results:collapsing results will remove duplicated guides
#'     found in overlapping ranges/regions.
#' @param remove_ambiguities Should spacer sequences with ambiguous 
#'     nucleotides ("N") be removed?
#'     TRUE by default. 
#' @param verbose Should messages be printed to console?
#'     TRUE by default.
#' 
#' @return A \linkS4class{GuideSet} object. 
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @details For input that is a character vector of genomic sequences,
#'     all identified guides are returned with their relative strand and
#'     coordinates. The input sequence is assumed to be on the positive
#'     strand starting at position 1. 
#' 
#' @examples
#' # Using custom sequence as input:
#' guides <- findSpacers("CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTAGCCGGGTGTTA")
#' 
#' # Exon-intro region of human KRAS specified 
#' # using a GRanges object:
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
#'     library(GenomicRanges)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' 
#'     gr_input <- GRanges(c("chr12"),
#'                         IRanges(start=25224014, end=25227007))
#'     genome(gr_input) <- "hg38"
#'     guideSet <- findSpacers(gr_input, bsgenome=bsgenome)
#' 
#'     # Designing guides for enCas12a nuclease:
#'     data(enAsCas12a, package="crisprBase")
#'     guideSet <- findSpacers(gr_input, 
#'                             bsgenome=bsgenome,
#'                             crisprNuclease=enAsCas12a)
#' }
#' 
#' @importFrom S4Vectors DataFrame mcols mcols<- 
#' @importFrom S4Vectors metadata metadata<- 
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics width strand strand<-
#' @importFrom BiocGenerics start end sort
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings matchPDict DNAStringSet reverseComplement
#' @importFrom IRanges resize trim Views flank 
#' @importFrom IRanges findOverlaps
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @importFrom GenomeInfoDb isCircular isCircular<- 
#' @importFrom GenomeInfoDb genome genome<- seqinfo seqinfo<-
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom crisprBase hasSpacerGap isCutting
#' @export
findSpacers <- function(x,
                        bsgenome=NULL,
                        crisprNuclease=NULL,
                        canonical=TRUE,
                        spacer_len=NULL,
                        both_strands=TRUE,
                        cut_offset=NULL,
                        strict_overlap=TRUE,
                        collapse=TRUE,
                        remove_ambiguities=TRUE,
                        verbose=FALSE
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    if (hasSpacerGap(crisprNuclease)){
        stop("CRISPR nucleases with spacer gaps are not ",
             "supported at the moment.")
    }
    if (!isCutting(crisprNuclease)){
        stop("CRISPR nucleases that are not cutting are not ",
             "supported at the moment.")
    }
    
    if (is(x, "GRanges")){
        gs <- .findSpacersFromGR(x,
                                 bsgenome=bsgenome,
                                 crisprNuclease=crisprNuclease,
                                 canonical=canonical,
                                 spacer_len=spacer_len,
                                 both_strands=both_strands,
                                 cut_offset=cut_offset,
                                 strict_overlap=strict_overlap,
                                 collapse=collapse,
                                 verbose=verbose)
    } else if (is(x, "character")){
        gs <- .findSpacersFromCustomSeq(x,
                                        crisprNuclease=crisprNuclease,
                                        canonical=canonical,
                                        spacer_len=spacer_len,
                                        both_strands=both_strands,
                                        cut_offset=cut_offset,
                                        collapse=collapse,
                                        verbose=verbose)
    } else {
        stop("x must be either a GRanges of a character vector. ")
    }
    gs <- BiocGenerics::sort(gs,
                             ignore.strand=TRUE)
    if (length(gs)>0){
        names(gs) <- paste0("spacer_", seq_along(gs))
    }
    if (remove_ambiguities){
        gs <- gs[!grepl("N", spacers(gs))]
    }
    return(gs)
}







#' @importFrom crisprBase spacerLength spacerLength<-
.findSpacersFromGR <- function(gr,
                               bsgenome=NULL,
                               crisprNuclease=NULL,
                               canonical=TRUE,
                               spacer_len=NULL,
                               both_strands=TRUE,
                               cut_offset=NULL,
                               strict_overlap=TRUE,
                               collapse=TRUE,
                               verbose=FALSE
){
    gr_original <- gr
    if (!is.null(spacer_len)){
        spacerLength(crisprNuclease) <- spacer_len
    } else {
        spacer_len <- spacerLength(crisprNuclease)
    }
    gr <- .validateGRanges(gr)
    gr <- .validateGRangesNames(gr)
    genome <- genome(gr)
    if (length(unique(genome))>1){
        stop("Multiple genomes were found for the input GRanges object.")
    }
    genome <- genome[1]
    if (is.null(bsgenome)){
        bsgenome  <- .getBSGenome(genome)
    } else {
        bsgenome <- .validateBSgenome(bsgenome)
        bsgenome_genome <- genome(bsgenome)[1]
        if (!is.na(genome)){
            if (bsgenome_genome != genome){
                stop("genome stored in the bsgenome object (",
                     bsgenome_genome, ") differs from genome provided ",
                     "in the input GRanges object (", genome, ").")
            }
        }
    }
    

    #-----------------------------------------------------------
    # Extracting DNA sequences into a DNAStringSet
    gr <- resize(gr,
                 fix="center",
                 width=width(gr)+2*spacer_len)
    gr <- trim(gr)
    strand(gr) <- "+"
    dna <- getSeq(bsgenome, gr)
    names(dna) <- names(gr)
    mcols(dna) <- as.data.frame(gr)[, c("seqnames", "start", "end")]
    metadata(dna)$genome  <- genome
    dna <- .validateInputFindSpacers(dna)
    #-----------------------------------------------------------


    gs <- .findSpacersFromDNAStringSet(dna,
                                       crisprNuclease=crisprNuclease,
                                       canonical=canonical,
                                       spacer_len=spacer_len,
                                       both_strands=both_strands,
                                       cut_offset=cut_offset,
                                       collapse=collapse,
                                       verbose=verbose)
    seqinfo(gs) <- seqinfo(gr_original)
    metadata(gs)[["bsgenome"]] <- bsgenome

    # Only keep sgRNAs with cut_site inside of the original GRanges:
    if (strict_overlap){
        gr.cutsite <- GRanges(seqnames(gs),
                              IRanges(start=gs$cut_site, width=1))
        hits <- findOverlaps(gr.cutsite,
                             gr_original,
                             ignore.strand=TRUE)
        hits <- sort(unique(queryHits(hits)))
        gs <- gs[hits]
    }
    return(gs)
}




#' @importFrom crisprBase spacerLength
.findSpacersFromCustomSeq <- function(x,
                                      crisprNuclease=NULL,
                                      canonical=TRUE,
                                      spacer_len=NULL,
                                      both_strands=TRUE,
                                      cut_offset=NULL,
                                      collapse=TRUE,
                                      verbose=FALSE
){

    if (!is.null(spacer_len)){
        spacerLength(crisprNuclease) <- spacer_len
    } else {
        spacer_len <- spacerLength(crisprNuclease)
    }

    # Transforming into a DNAStringSet
    x <- .validateCustomSeqNames(x)
    x <- .validateDNACharacterVariable(seq=x,
                                       argument="x",
                                       nullOk=FALSE,
                                       exactBases=FALSE)
    dna <- DNAStringSet(x)
    metadata(dna)$genome <- "custom"
    mcols(dna) <- data.frame(seqnames=names(dna),
                             start=1,
                             end=width(dna))
    dna <- .validateInputFindSpacers(dna)
    

    gs <- .findSpacersFromDNAStringSet(dna,
                                       crisprNuclease=crisprNuclease,
                                       canonical=canonical,
                                       spacer_len=spacer_len,
                                       both_strands=both_strands,
                                       cut_offset=cut_offset,
                                       collapse=collapse,
                                       verbose=verbose)
    metadata(gs)[["custom_seq"]] <- dna

    #Taking care of seqinfo:
    lens <- width(metadata(gs)[["custom_seq"]])
    names(lens) <- names(metadata(gs)[["custom_seq"]])
    if (length(gs)>0){
        seqlengths(gs) <- lens
        isCircular(gs) <- rep(FALSE, length(lens))
    }
    genome(gs) <- "custom"
    return(gs)
}







#' @importFrom crisprBase pams spacerSide
.findSpacersFromDNAStringSet <- function(dna,
                                         crisprNuclease=NULL,
                                         canonical=TRUE,
                                         spacer_len=NULL,
                                         both_strands=TRUE,
                                         cut_offset=NULL,
                                         collapse=TRUE,
                                         verbose=FALSE
){

    pams <- crisprBase::pams(crisprNuclease,
                             primary=canonical)
    pams <- DNAStringSet(pams)

        
    getPerSequenceHits <- function(pams,
                                   dna_string,
                                   strand="+",
                                   name){
        pam_ir <- matchPDict(pams,dna_string) 
        pam_ir <- Reduce(c, pam_ir)
        if (spacerSide(crisprNuclease)=="5prime"){
            spacer_ir <- flank(pam_ir,
                               width=spacer_len,
                               start=TRUE)
        } else {
            spacer_ir <- flank(pam_ir,
                               width=spacer_len,
                               start=FALSE)
        }
        valid_start <- BiocGenerics::start(spacer_ir)>0
        valid_end   <- BiocGenerics::end(spacer_ir)<=length(dna_string)
        valid <- valid_start & valid_end
        pam_ir    <- pam_ir[valid]
        spacer_ir <- spacer_ir[valid]
        sequences_pam    <- DNAStringSet(Views(dna_string, pam_ir))
        sequences_spacer <- DNAStringSet(Views(dna_string, spacer_ir))
        list(pam=sequences_pam,
             spacer=sequences_spacer,
             pam_site=BiocGenerics::start(pam_ir),
             region=rep(name, length(pam_ir)),
             strand=rep(strand, length(pam_ir)))
    }

    # Getting hits:
    hits <- lapply(seq_along(dna), function(i){
        getPerSequenceHits(pams,
                           dna[[i]],
                           strand="+",
                           name=names(dna)[i])
    })
    if (both_strands){
        dna_rev  <- reverseComplement(dna)
        hits_rev <- lapply(seq_along(dna_rev), function(i){
            getPerSequenceHits(pams,
                               dna_rev[[i]],
                               strand="-",
                               name=names(dna)[i])
        })
        hits <- c(hits,hits_rev)
    }

    # Merging everything together
    spacer <- lapply(hits, function(x) x$spacer)
    out    <- DataFrame(spacer=do.call(c, spacer))
    out$pam      <- do.call(c, lapply(hits, function(x) x$pam))
    out$pam_site <- do.call(c, lapply(hits, function(x) x$pam_site))
    out$strand   <- do.call(c, lapply(hits, function(x) x$strand))
    out$region   <- do.call(c, lapply(hits, function(x) x$region))

    # Converting coordinate to fwd strand
    # for gRNAs located on rev strand:
    wh <- which(out$strand=="-")
    if (length(wh)>0){
        lens=width(dna)[match(out$region, names(dna))]
        out$pam_site[wh] <- lens[wh] - out$pam_site[wh] + 1
    }
    
    # Transforming coordinates to original coordinates system:
    wh <- match(out$region, names(dna))
    out$pam_site <- out$pam_site + mcols(dna)$start[wh] -1 
    out$chr <- mcols(dna)$seqnames[wh]

    # Creating GuideSet
    gs <- GuideSet(spacers=out$spacer,
                   pams=out$pam,
                   seqnames=out$chr,
                   pam_site=out$pam_site,
                   strand=out$strand,
                   CrisprNuclease=crisprNuclease,
                   genome=metadata(dna)$genome)
    mcols(gs)$region <- out$region
    if (length(gs)>0){
        names(gs) <- paste0("spacer_", seq_along(gs))
    }
    if ("custom_seq" %in% names(metadata(dna))){
        metadata(gs)[["custom_seq"]] <- metadata(dna)[["custom_seq"]]
    }
    
    if (metadata(gs)$genome=="custom"){
        collapse=FALSE
    }
    if (collapse){
        mcols(gs)[["region"]] <- NULL
        df <- as.data.frame(gs)[, c("seqnames", "start", "strand")]
        gs <- gs[!duplicated(df)]
        if (length(gs)>0){
            names(gs) <- paste0("spacer_", seq_along(gs))
        }
    }
    cut_site <- getCutSiteFromPamSite(pam_site=pamSites(gs),
                                      strand=as.character(strand(gs)),
                                      crisprNuclease=crisprNuclease(gs),
                                      cut_offset=cut_offset)
    mcols(gs)[["cut_site"]] <- cut_site
    return(gs)
}


# .convertGRToGuideSet <- function(gr){
#     guideSet <- GuideSet(CrisprNuclease=metadata(gr)$nuclease,
#                          genome=metadata(gr)$genome,
#                          pam_site=mcols(gr)$pam_site,
#                          strand=as.character(strand(gr)),
#                          pams=as.character(mcols(gr)$pam),
#                          spacers=as.character(mcols(gr)$spacer),
#                          seqnames=as.character(seqnames(gr)))
#     mcols(guideSet)$cut_site <- mcols(gr)$cut_site
#     names(guideSet) <- names(gr)
#     gr_cols <- colnames(mcols(gr))
#     gs_cols <- colnames(mcols(guideSet))
#     other_cols <- setdiff(gr_cols, gs_cols)
#     if (length(other_cols)>0){
#         mcols(guideSet)[,other_cols] <- mcols(gr)[,other_cols]
#     }
#     seqinfo(guideSet) <- seqinfo(gr)
#     return(guideSet)
# }






