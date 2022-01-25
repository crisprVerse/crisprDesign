#' @title Find CRISPR gRNA spacer sequences in a set of DNA sequences.
#' @description Returns all possible, valid gRNA sequences from either a
#'    \linkS4class{GRanges} object or a character vector of genomic sequences,
#'    for a given CRISPR nuclease.
#' 
#' @param x Either a \linkS4class{GRanges}, a \linkS4class{DNAStringSet}, or a
#'     \linkS4class{DNAString} object, or a character vector of genomic
#'     sequences. See details.
#' @param bsgenome \linkS4class{BSgenome} object from which to extract
#'     sequences if a \linkS4class{GRanges} object is provided as input. 
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param canonical Whether to return only guide sequences having canonical
#'     PAM sequences. If TRUE (default), only PAM sequences with the highest
#'     weights stored in the \code{crisprNuclease} object will be considered.
#' @param spacer_len Length of spacer to be used if different from
#'     default length specified in \code{crisprNuclease}.
#' @param both_strands Should both strands be considered to search for 
#'     protospacer sequences? TRUE by default.
#' @param cut_offset Distance in nucleotides between \code{pam_site}
#'     and \code{cut_site} if different from default offset specified by
#'     \code{crisprNuclease}.
#' @param strict_overlap If \code{TRUE} (default), only include guides
#'     having the cut site in the input range. Ignored when \code{x} is not a
#'     \linkS4class{GRanges} object.
#' @param remove_ambiguities Should spacer sequences containing ambiguous
#'     nucleotides (not "A", "C", "G", or "T") be removed? TRUE by default.
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
#'     # Designing guides for enAsCas12a nuclease:
#'     data(enAsCas12a, package="crisprBase")
#'     guideSet <- findSpacers(gr_input, 
#'                             canonical=FALSE,
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
#' @importFrom crisprBase spacerLength<-
#' @export
findSpacers <- function(x,
                        crisprNuclease=NULL,
                        bsgenome=NULL,
                        canonical=TRUE,
                        both_strands=TRUE,
                        spacer_len=NULL,
                        cut_offset=NULL,
                        strict_overlap=TRUE,
                        remove_ambiguities=TRUE
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    .checkCrisprNucleaseForSupportedFeatures(crisprNuclease)
    for (i in (c("canonical", "both_strands", "remove_ambiguities"))){
        .checkSingleBoolean(i, get(i))
    }
    .checkSingleInteger("spacer_len", spacer_len, sign="positive")
    if (!is.null(spacer_len)){
        crisprBase::spacerLength(crisprNuclease) <- spacer_len
    }
    .checkSingleInteger("cut_offset", cut_offset)
    
    dna <- .asDNAStringSet(x,
                           bsgenome=bsgenome,
                           crisprNuclease=crisprNuclease,
                           both_strands=both_strands)
    gs <- .findSpacersFromDNAStringSet(dna=dna,
                                       crisprNuclease=crisprNuclease,
                                       canonical=canonical,
                                       both_strands=both_strands,
                                       cut_offset=cut_offset,
                                       dnaFromGenome=.isGRanges(x))
    gs <- .addMetadata(gs=gs,
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






# update when features are supported
#' @importFrom crisprBase hasSpacerGap isRnase isCutting
.checkCrisprNucleaseForSupportedFeatures <- function(crisprNuclease
){
    if (crisprBase::hasSpacerGap(crisprNuclease)){
        stop("CRISPR nucleases with spacer gaps are not ",
             "supported at the moment.")
    }
    if (crisprBase::isRnase(crisprNuclease)){
        stop("RNA-targeting CRISPR nucleases are not ",
             "supported at the moment.")
    }
    if (!crisprBase::isCutting(crisprNuclease)){
        stop("CRISPR nucleases that are not cutting are not ",
             "supported at the moment.")
    }
    invisible(NULL)
}



#' @importFrom methods is
#' @importFrom BiocGenerics strand strand<-
.asDNAStringSet <- function(x,
                            bsgenome,
                            crisprNuclease,
                            both_strands
){
    if (.isGRanges(x)){
        if (any(as.character(BiocGenerics::strand(x)) == "*")){
            ambiguousStrand <-as.character(BiocGenerics::strand(x)) == "*"
            BiocGenerics::strand(x)[ambiguousStrand] <- "+"
            if (!both_strands){
                revRanges <- x[ambiguousStrand]
                BiocGenerics::strand(revRanges) <- "-"
                x <- c(x, revRanges)
            }
        }
        x <- .GRanges2DNAStringSet(x,
                                   bsgenome=bsgenome,
                                   crisprNuclease=crisprNuclease)
    } else if (methods::is(x, "DNAString") || is.vector(x, mode="character")){
        x <- .string2DNAStringSet(x)
    }
    if (!methods::is(x, "DNAStringSet")){
        stop("Value type for 'x' not recognized; see ?findSpacers")
    }
    return(x)
}






#' @importFrom GenomeInfoDb genome seqnames
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges resize trim
#' @importFrom BiocGenerics width start end strand
#' @importFrom S4Vectors mcols<- metadata<- DataFrame
#' @importFrom crisprBase spacerLength
.GRanges2DNAStringSet <- function(x,
                                  bsgenome,
                                  crisprNuclease
){
    x <- .validateGRanges(x)
    x <- .validateGRangesNames(x)
    genome <- unique(GenomeInfoDb::genome(x))
    if (length(genome) > 1){
        stop("Multiple genomes found for the input GRanges object.")
    }
    if (is.null(bsgenome)){
        bsgenome  <- .getBSGenome(genome)
    } else {
        bsgenome <- .validateBSgenome(bsgenome)
        bsgenome_genome <- unique(GenomeInfoDb::genome(bsgenome))
        if (!is.na(genome) && genome != bsgenome_genome){
            stop("genome stored in the bsgenome object (",
                 bsgenome_genome, ") differs from genome provided ",
                 "in the input GRanges object (", genome, ").")
        }
    }
    newWidth <- BiocGenerics::width(x) +
        2 * crisprBase::spacerLength(crisprNuclease)
    x <- GenomicRanges::resize(x,
                               fix="center",
                               width=newWidth)
    x <- GenomicRanges::trim(x)
    dna <- BSgenome::getSeq(bsgenome, x)
    mcols <- S4Vectors::DataFrame(
        seqnames=as.character(GenomeInfoDb::seqnames(x)),
        start=BiocGenerics::start(x),
        end=BiocGenerics::end(x),
        strand=as.character(BiocGenerics::strand(x)))
    S4Vectors::mcols(dna) <- mcols
    
    S4Vectors::metadata(dna)$genome  <- genome
    return(dna)
}




#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors DataFrame mcols<- metadata<-
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb genome
.string2DNAStringSet <- function(x
){
    dna <- Biostrings::DNAStringSet(x)
    regionNames <- names(dna)
    if (is.null(regionNames)){
        regionNames <- paste0("region_", length(dna), recycle0=TRUE)
    }
    names(dna) <- regionNames
    
    mcols <- S4Vectors::DataFrame(seqnames=regionNames,
                                  start=1,
                                  end=BiocGenerics::width(dna),
                                  strand="+")
    S4Vectors::mcols(dna) <- mcols
    S4Vectors::metadata(dna)$genome <- "custom"
    return(dna)
}




#' @importFrom crisprBase pams
#' @importFrom Biostrings reverseComplement
#' @importFrom S4Vectors mcols mcols<- bindROWS
#' @importFrom BiocGenerics width
.findSpacersFromDNAStringSet <- function(dna,
                                         crisprNuclease,
                                         canonical,
                                         both_strands,
                                         cut_offset,
                                         dnaFromGenome
){
    pams <- crisprBase::pams(crisprNuclease, primary=canonical)
    if (both_strands){
        # more efficient to "reverseComplement" crisprNuclease for large seqs?
        revComp <- Biostrings::reverseComplement(dna)
        revCompStrand <- chartr("+-", "-+", S4Vectors::mcols(revComp)$strand)
        S4Vectors::mcols(revComp)$strand <- revCompStrand
        S4Vectors::mcols(dna)$revComp <- FALSE
        S4Vectors::mcols(revComp)$revComp <- TRUE
        dna <- c(dna, revComp)
    } else {
        S4Vectors::mcols(dna)$revComp <- FALSE
    }
    hits <- lapply(seq_along(dna), function(i){
        .spacersPerSequence(seq=dna[i],
                            pams=pams,
                            crisprNuclease=crisprNuclease)
    })
    hits <- Reduce(S4Vectors::bindROWS, hits)
    if (both_strands){
        # fix relative pam_site for complement strand
        regions <- hits$region[hits$revComp]
        seqWidths <- BiocGenerics::width(dna[regions])
        newPamSites <- seqWidths - hits$pam_site[hits$revComp] + 1
        hits$pam_site[hits$revComp] <- newPamSites
        hits$revComp <- NULL
    }
    
    if (dnaFromGenome){
        # flip pam_sites on negative strand
        revStrand <- S4Vectors::mcols(dna)$strand == "-" &
            !S4Vectors::mcols(dna)$revComp
        revRegions <- row.names(mcols(dna))[revStrand]
        revSpacers <- hits$region == revRegions
        regionWidths <- hits$region[revSpacers]
        regionWidths <- width(dna[regionWidths])
        newPamSites <- regionWidths - hits$pam_site[revSpacers] + 1
        hits$pam_site[revSpacers] <- newPamSites
        # transform pam_sites to genomic coordinates
        newStarts <- S4Vectors::mcols(dna)[hits$region, ]
        newStarts <- newStarts$start
        hits$pam_site <- newStarts + hits$pam_site - 1
    }
    
    indices <- match(hits$region, names(dna))
    hits$chr <- S4Vectors::mcols(dna)$seqnames[indices]
    
    S4Vectors::mcols(dna)$revComp <- NULL
    
    gs <- GuideSet(protospacers=hits$spacer,
                   pams=hits$pam,
                   seqnames=hits$chr,
                   pam_site=hits$pam_site,
                   strand=hits$strand,
                   CrisprNuclease=crisprNuclease,
                   genome=metadata(dna)$genome)
    genome(gs) <- metadata(dna)$genome # incorporate into GuideSet constructor
    cut_site <- getCutSiteFromPamSite(pam_site=pamSites(gs),
                                      strand=as.character(strand(gs)),
                                      crisprNuclease=crisprNuclease(gs),
                                      cut_offset=cut_offset)
    mcols(gs)[["cut_site"]] <- cut_site
    mcols(gs)[['region']] <- hits$region
    metadata(gs)[["CrisprNuclease"]] <- crisprNuclease
    return(gs)
}




#' @importFrom Biostrings matchPDict extractAt
#' @importFrom crisprBase spacerLength pamSide
#' @importFrom GenomicRanges flank
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
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
    DataFrame(pam=pamSequences,
              spacer=spacerSequences,
              pam_site=BiocGenerics::start(pamRanges),
              region=rep(names(seq), validCount),
              strand=rep(S4Vectors::mcols(seq)$strand, validCount),
              revComp=rep(S4Vectors::mcols(seq)$revComp, validCount))
}




#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqinfo seqinfo<-
#' @importFrom GenomeInfoDb genome
#' @importFrom S4Vectors metadata<-
.addMetadata <- function(gs,
                         x,
                         dna,
                         bsgenome
                         
){
    if (.isGRanges(x)){
        GenomeInfoDb::seqlevels(gs) <- GenomeInfoDb::seqlevels(x)
        GenomeInfoDb::seqinfo(gs) <- GenomeInfoDb::seqinfo(x)
        genome <- GenomeInfoDb::genome(GenomeInfoDb::seqinfo(gs))
        S4Vectors::metadata(gs)[["genome"]] <- unique(genome)
        if (is.null(bsgenome)){
            genome <- unique(GenomeInfoDb::genome(x))
            bsgenome  <- .getBSGenome(genome)
        }
        S4Vectors::metadata(gs)[["bsgenome"]] <- bsgenome
    } else {
        customSeqInfo <- Seqinfo(seqnames=names(dna),
                                 seqlengths=width(dna),
                                 isCircular=NA,
                                 genome="custom")
        GenomeInfoDb::seqinfo(gs) <- customSeqInfo
        S4Vectors::metadata(gs)[["genome"]] <- "custom"
        S4Vectors::metadata(gs)[["custom_seq"]] <- dna
    }
    return(gs)
}



#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
.applyStrictOverlap <- function(gs,
                                x,
                                strict_overlap
){
    .checkSingleBoolean("strict_overlap", strict_overlap)
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
        if (crisprBase::isDnase(crisprNuclease)){
            bases <- Biostrings::DNA_BASES
        } else if (crisprBase::isRnase(crisprNuclease)){
            bases <- Biostrings::RNA_BASES
        }
        bases <- paste0(bases, collapse="")
        pattern <- paste0("^[", bases, "]+$")
        guideSet <- guideSet[grepl(pattern, spacers(guideSet))]
    }
    return(guideSet)
}
