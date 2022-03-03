#' @title Find CRISPR gRNA spacer sequences from a set of DNA sequences.
#' @description Returns all possible, valid gRNA sequences for a given CRISPR
#'    nuclease from either a \linkS4class{GRanges} object or a set of
#'    sequence(s) contained in either a \linkS4class{DNAStringSet},
#'    \linkS4class{DNAString} or character vector of genomic sequences.
#' 
#' @param x Either a \linkS4class{GRanges}, a \linkS4class{DNAStringSet}, or a
#'     \linkS4class{DNAString} object, or a character vector of genomic
#'     sequences. See details.
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
#' @return A \linkS4class{GuideSet} object. 
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @details If \code{x} is a \linkS4class{GRanges} object then a
#'     \linkS4class{BSgenome} must be supplied to \code{bsgenome}, from which
#'     the genomic sequence is obtained, unless the \code{bsgenome} can be
#'     inferred from \code{genome(x)}, for example, \code{"hg38"}. Otherwise,
#'     all supplied sequences are treated as the \code{"+"} strands of
#'     chromosomes in a \code{"custom"} genome.
#'     
#'     Ranges or sequences in \code{x} may contain names where permitted. These
#'     names are stored in \code{region} in the \code{mcols} of the output,
#'     and as \code{seqnames} of the output if \code{x} is not a
#'     \linkS4class{GRanges} object. If not \code{NULL}, \code{names(x)} must
#'     be unique, otherwise ranges or sequences are enumerated with the
#'     \code{"region_"} prefix.
#'     
#'     When \code{x} is a \linkS4class{GRanges}, the \code{*} strand is
#'     interpreted as both strands. Consequently, the \code{both_strands}
#'     argument has no effect on such ranges.
#'     
#' 
#' @examples
#' # Using custom sequence as input:
#' my_seq <- c(my_seq="CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTAGCCGGGTGTTA")
#' guides <- findSpacers(my_seq)
#' 
#' # Exon-intro region of human KRAS specified 
#' # using a GRanges object:
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
#'     library(GenomicRanges)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' 
#'     gr_input <- GRanges(c("chr12"),
#'                         IRanges(start=25224014, end=25227007))
#'     guideSet <- findSpacers(gr_input, bsgenome=bsgenome)
#' 
#'     # Designing guides for enAsCas12a nuclease:
#'     data(enAsCas12a, package="crisprBase")
#'     guideSet <- findSpacers(gr_input, 
#'                             canonical=FALSE,
#'                             bsgenome=bsgenome,
#'                             crisprNuclease=enAsCas12a)
#' }
#' 
#' 
#' @importFrom BiocGenerics sort
#' @export
findSpacers <- function(x,
                        crisprNuclease=NULL,
                        bsgenome=NULL,
                        canonical=TRUE,
                        both_strands=TRUE,
                        spacer_len=NULL,
                        strict_overlap=TRUE,
                        remove_ambiguities=TRUE
){
    if (.isGRanges(x)){
        targetOrigin <- "bsgenome"
        customSequences <- NULL
    } else {
        targetOrigin <- "customSequences"
        customSequences <- x
    }
    crisprNuclease <- .setupCrisprNucleaseForFindSpacers(crisprNuclease=crisprNuclease,
                                                         spacer_len=spacer_len)
    for (i in (c("canonical", "both_strands", "remove_ambiguities"))){
        .checkBoolean(i, get(i))
    }
    if (isRnase(crisprNuclease)){
        both_strands=FALSE
    }
    dna <- .asDNAStringSet(x,
                           bsgenome=bsgenome,
                           crisprNuclease=crisprNuclease,
                           both_strands=both_strands)
    gs <- .findSpacersFromDNAStringSet(dna=dna,
                                       bsgenome=bsgenome,
                                       customSequences=customSequences,
                                       crisprNuclease=crisprNuclease,
                                       canonical=canonical,
                                       targetOrigin=targetOrigin)
    gs <- .cleanSeqInfo(gs=gs,
                        x=x,
                        dna=dna,
                        bsgenome=bsgenome)
    gs <- .applyStrictOverlap(gs=gs,
                              x=x,
                              strict_overlap=strict_overlap)
    gs <- .removeAmbiguities(guideSet=gs,
                             crisprNuclease=crisprNuclease,
                             remove_ambiguities=remove_ambiguities)
    gs <- BiocGenerics::sort(gs, ignore.strand=TRUE)
    names(gs) <- paste0("spacer_", seq_along(gs), recycle0=TRUE)
    return(gs)
}



#' @importFrom crisprBase spacerLength<-
.setupCrisprNucleaseForFindSpacers <- function(crisprNuclease,
                                               spacer_len
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    .checkCrisprNucleaseForSupportedFeatures(crisprNuclease)
    .checkSingleInteger("spacer_len", spacer_len, sign="positive")
    if (!is.null(spacer_len)){
        crisprBase::spacerLength(crisprNuclease) <- spacer_len
    }
    return(crisprNuclease)
}


# update when features are supported
#' @importFrom crisprBase hasSpacerGap isRnase isCutting
.checkCrisprNucleaseForSupportedFeatures <- function(crisprNuclease
){
    if (crisprBase::hasSpacerGap(crisprNuclease)){
        stop("CRISPR nucleases with spacer gaps are not ",
             "supported at the moment.")
    }
    #if (crisprBase::isRnase(crisprNuclease)){
        #stop("RNA-targeting CRISPR nucleases are not ",
    #         "supported at the moment.")
    #}
    #if (!crisprBase::isCutting(crisprNuclease)){
    #    stop("CRISPR nucleases that are not cutting are not ",
    #         "supported at the moment.")
    #}
    invisible(NULL)
}


#' @importFrom methods is
.asDNAStringSet <- function(x,
                            bsgenome,
                            crisprNuclease,
                            both_strands
){
    if (!is.null(names(x)) && any(duplicated(names(x)))){
        stop("Provided names for 'x' must be unique")
    }
    if (any(is.na(names(x)))){
        stop("Names for 'x' cannot be NA")
    }
    if (.isGRanges(x)){
        x <- .GRanges2DNAStringSet(x,
                                   bsgenome=bsgenome,
                                   both_strands=both_strands,
                                   crisprNuclease=crisprNuclease)
    } else if (methods::is(x, "DNAString") || 
               methods::is(x, "DNAStringSet") ||
               is.vector(x, mode="character")){
        x <- .string2DNAStringSet(x,
                                  both_strands=both_strands)
    } else {
        stop("Value type for 'x' not recognized; see ?findSpacers")
    }
    return(x)
}


#' @importFrom GenomeInfoDb genome seqlevels seqlevels<- seqinfo seqinfo<-
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BSgenome getSeq
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors DataFrame mcols<- metadata<-
.GRanges2DNAStringSet <- function(x,
                                  bsgenome,
                                  both_strands,
                                  crisprNuclease
){
    genome <- unique(GenomeInfoDb::genome(x))
    if (length(genome) > 1){
        stop("Multiple genomes found for the input GRanges object.")
    }
    bsgenome <- .bsgenome4GrangesInput(bsgenome=bsgenome,
                                       genome=genome)
    GenomeInfoDb::seqlevels(x) <- GenomeInfoDb::seqlevels(bsgenome)
    GenomeInfoDb::seqinfo(x) <- GenomeInfoDb::seqinfo(bsgenome)
    x <- .assignRegionNames(x)
    x <- .resolveRegionStrands(x, both_strands=both_strands)
    x <- .expandGrangesBySpacerLength(x, crisprNuclease=crisprNuclease)
    
    dna <- BSgenome::getSeq(bsgenome, x)
    S4Vectors::mcols(dna) <- S4Vectors::DataFrame(
        seqnames=as.character(GenomeInfoDb::seqnames(x)),
        start=BiocGenerics::start(x),
        end=BiocGenerics::end(x),
        strand=as.character(BiocGenerics::strand(x)))
    return(dna)
}


#' @importFrom GenomeInfoDb genome
.bsgenome4GrangesInput <- function(bsgenome,
                                   genome
){
    if (is.null(bsgenome)){
        stop("bsgenome must be provided.")
    } else {
        .isBSGenome(bsgenome)
        bsgenome_genome <- unique(GenomeInfoDb::genome(bsgenome))
        if (!is.na(genome) && genome != bsgenome_genome){
            stop("genome stored in the bsgenome object (",
                 bsgenome_genome, ") differs from genome provided ",
                 "in the input GRanges object (", genome, ").")
        }
    }
    return(bsgenome)
}


.assignRegionNames <- function(x){
    regionNames <- names(x)
    if (is.null(regionNames) ||
        all(is.na(regionNames)) ||
        all(regionNames == "")){
        regionNames <- paste0("region_", seq_along(x), recycle0=TRUE)
    }
    names(x) <- regionNames
    return(x)
}


#' @importFrom BiocGenerics strand strand<- invertStrand
.resolveRegionStrands <- function(x,
                                  both_strands
){
    if (any(BiocGenerics::strand(x) == "*")){
        ambiguousStrand <- as.character(BiocGenerics::strand(x)) == "*"
        BiocGenerics::strand(x)[ambiguousStrand] <- "+"
        if (!both_strands){
            revStrand <- BiocGenerics::invertStrand(x[ambiguousStrand])
            x <- c(x, revStrand)
        }
    }
    if (both_strands){
        revStrand <- BiocGenerics::invertStrand(x)
        x <- c(x, revStrand)
    }
    return(x)
}


#' @importFrom crisprBase spacerLength
#' @importFrom BiocGenerics width
#' @importFrom GenomicRanges resize trim
.expandGrangesBySpacerLength <- function(x,
                                         crisprNuclease
){
    spacer_length <- crisprBase::spacerLength(crisprNuclease)
    new_width <- BiocGenerics::width(x) + 2 * spacer_length
    x <- GenomicRanges::resize(x,
                               fix="center",
                               width=new_width)
    x <- GenomicRanges::trim(x)
    return(x)
}


#' @importFrom Biostrings DNAStringSet reverseComplement
#' @importFrom S4Vectors mcols mcols<- DataFrame metadata<-
#' @importFrom BiocGenerics width invertStrand
.string2DNAStringSet <- function(x,
                                 both_strands
){
    dna <- Biostrings::DNAStringSet(x)
    dna <- .assignRegionNames(dna)
    S4Vectors::mcols(dna) <- S4Vectors::DataFrame(
        seqnames=names(dna),
        start=1,
        end=BiocGenerics::width(dna),
        strand="+")
    S4Vectors::metadata(dna)$genome <- "custom"
    if (both_strands){
        revComp <- Biostrings::reverseComplement(dna)
        mcols <- S4Vectors::mcols(revComp)
        S4Vectors::mcols(revComp) <- BiocGenerics::invertStrand(mcols)
        dna <- c(dna, revComp)
    }
    return(dna)
}


#' @importFrom crisprBase pams
#' @importFrom S4Vectors mcols<- metadata metadata<- bindROWS
#' @importFrom GenomeInfoDb genome<-
#' @importFrom BiocGenerics strand
.findSpacersFromDNAStringSet <- function(dna,
                                         bsgenome,
                                         customSequences,
                                         crisprNuclease,
                                         canonical,
                                         targetOrigin
){

    pams <- crisprBase::pams(crisprNuclease, primary=canonical)
    hits <- lapply(seq_along(dna), function(i){
        .spacersPerSequence(seq=dna[i],
                            pams=pams,
                            crisprNuclease=crisprNuclease)
    })
    hits <- Reduce(S4Vectors::bindROWS, hits)
    hits <- .applyHitsCoordinates(hits=hits,
                                  dna=dna)
    
    gs <- GuideSet(protospacers=hits$spacer,
                   pams=hits$pam,
                   seqnames=hits$chr,
                   pam_site=hits$pam_site,
                   strand=hits$strand,
                   CrisprNuclease=crisprNuclease,
                   targetOrigin=targetOrigin,
                   bsgenome=bsgenome,
                   customSequences=customSequences)
    cut_site <- getCutSiteFromPamSite(pam_site=pamSites(gs),
                                      strand=as.character(BiocGenerics::strand(gs)),
                                      crisprNuclease=crisprNuclease(gs))
    S4Vectors::mcols(gs)[["cut_site"]] <- cut_site
    S4Vectors::mcols(gs)[['region']] <- hits$region
    S4Vectors::metadata(gs)[["CrisprNuclease"]] <- crisprNuclease
    return(gs)
}




#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqinfo seqinfo<- genome
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importFrom S4Vectors metadata<-
#' @importFrom BiocGenerics width
.cleanSeqInfo <- function(gs,
                          x,
                          dna,
                          bsgenome
){
    if (.isGRanges(x)){
        GenomeInfoDb::seqlevels(gs) <- GenomeInfoDb::seqlevels(bsgenome)
        GenomeInfoDb::seqinfo(gs) <- GenomeInfoDb::seqinfo(bsgenome)
    } else {
        dna_names <- unique(names(dna))
        dna_lengths <- BiocGenerics::width(dna[dna_names])
        customSeqInfo <- GenomeInfoDb::Seqinfo(seqnames=dna_names,
                                               seqlengths=dna_lengths,
                                               isCircular=NA,
                                               genome="custom")
        GenomeInfoDb::seqinfo(gs) <- customSeqInfo
    }
    return(gs)
}




#' @importFrom Biostrings matchPDict extractAt
#' @importFrom crisprBase spacerLength pamSide
#' @importFrom GenomicRanges flank
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors DataFrame mcols
.spacersPerSequence <- function(seq,
                                pams,
                                crisprNuclease
){
    seqString <- unlist(seq)
    pamRanges <- Biostrings::matchPDict(pams, seqString) 
    pamRanges <- Reduce(c, pamRanges)
    flankWidth <- crisprBase::spacerLength(crisprNuclease)
    flankAtStart <- crisprBase::pamSide(crisprNuclease) == "3prime"
    spacerRanges <- GenomicRanges::flank(pamRanges,
                                         width=flankWidth,
                                         start=flankAtStart)
    validStart <- BiocGenerics::start(spacerRanges) > 0
    validEnd <- BiocGenerics::end(spacerRanges) <= length(seqString)
    valid <- validStart & validEnd
    pamRanges    <- pamRanges[valid]
    spacerRanges <- spacerRanges[valid]
    pamSequences <- Biostrings::extractAt(seqString, pamRanges)
    spacerSequences <- Biostrings::extractAt(seqString, spacerRanges)
    validCount <- sum(valid)
    spacers <- S4Vectors::DataFrame(
        pam=pamSequences,
        spacer=spacerSequences,
        pam_site=BiocGenerics::start(pamRanges),
        region=rep(names(seq), validCount),
        strand=rep(S4Vectors::mcols(seq)$strand, validCount))
    return(spacers)
}


#' @importFrom BiocGenerics width
#' @importFrom S4Vectors mcols
.applyHitsCoordinates <- function(hits,
                                  dna
){
    revStrand <- hits$strand == "-"
    regionLengths <- BiocGenerics::width(dna[hits$region[revStrand]])
    new_pam_site <- regionLengths - hits$pam_site[revStrand] + 1
    hits$pam_site[revStrand] <- new_pam_site
    
    indices <- match(hits$region, names(dna))
    regionStarts <- S4Vectors::mcols(dna)$start[indices]
    hits$pam_site <- regionStarts + hits$pam_site - 1
    hits$chr <- S4Vectors::mcols(dna)$seqnames[indices]
    return(hits)
}





#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
.applyStrictOverlap <- function(gs,
                                x,
                                strict_overlap
){
    .checkBoolean("strict_overlap", strict_overlap)
    if (.isGRanges(x) && strict_overlap){
        cut_sites <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(gs),
                                            IRanges::IRanges(start=gs$cut_site,
                                                             width=1))
        hits <- GenomicRanges::findOverlaps(cut_sites,
                                            x,
                                            ignore.strand=TRUE)
        hits <- unique(S4Vectors::queryHits(hits))
        gs <- gs[hits]
    }
    return(gs)
}


#' @importFrom crisprBase isDnase isRnase
#' @importFrom Biostrings DNA_BASES RNA_BASES
.removeAmbiguities <- function(guideSet,
                               crisprNuclease,
                               remove_ambiguities
){
    if (remove_ambiguities){
        #if (crisprBase::isDnase(crisprNuclease)){
        #    bases <- Biostrings::DNA_BASES
        #} else if (crisprBase::isRnase(crisprNuclease)){
        #    bases <- Biostrings::RNA_BASES
        #}
        bases <- Biostrings::DNA_BASES
        bases <- paste0(bases, collapse="")
        pattern <- paste0("^[", bases, "]+$")
        guideSet <- guideSet[grepl(pattern, spacers(guideSet))]
    }
    return(guideSet)
}

