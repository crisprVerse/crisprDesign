#' @title Functions for finding and characterizing on- and off-targets of
#'     spacer sequences.
#' 
#' @description Functions for finding and characterizing on- and off-targets of
#'     spacer sequences.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param aligner Which genomic alignment method should be used?
#'     Must be one of "bowtie", "bwa", and "biostrings".
#'    "bowtie" by default. Note that "bwa" is not availble for 
#'     Windows machines. 
#' @param columnName String specifying the columm name storing the alignments
#'     in \code{mcols(guideSet)}. "alignments" by default.
#' @param addSummary Should summary columns be added to \code{guideSet}?
#'     TRUE by default.
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}} for annotating
#'     on-target and off-target alignments using gene annotation.
#' @param tssObject A \linkS4class{GRanges} object specifying TSS coordinates.
#' @param spacers Character vector of gRNA spacer sequences.
#'     All sequences must be equal in length.
#' @param custom_seq Optional string specifying the target DNA sequence for
#'     the search space. This will limit the off-target
#'     search to the specified custom sequence. 
#' @param aligner_index String specifying bowtie or BWA index.
#'     Must be provided when \code{aligner} is either \code{"bowtie"} or
#'     \code{"bwa"}.
#' @param bsgenome A \linkS4class{BSgenome} object from which to extract
#'     sequences if a \linkS4class{GRanges} object is provided as input. 
#' @param n_mismatches Maximum number of mismatches permitted between guide RNA
#'     and genomic DNA.
#' @param n_max_alignments Maximum number of alignments to report by bowtie 
#'     for each spacer. Effectively set to \code{Inf} when \code{allPossible}
#'     is \code{TRUE}.
#' @param all_alignments Should all all possible alignments be returned?
#'     FALSE by defaule.
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param canonical Should only alignments corresponding to canonical
#'     PAM sequences be returned? TRUE by default.
#' @param ignore_pam If TRUE, will return all matches regardless of
#'     PAM sequence. FALSE by default. 
#' @param standard_chr_only Should only standard chromosomes be considered?
#'     TRUE by default.
#' @param tss_window Window size of promoters upstream of gene TSS to search
#'     for overlap with spacer sequence. Must be a numeric vector of length 2:
#'     upstream limit and downstream limit. Default is \code{c(-500, 500)},
#'     which includes 500bp upstream and downstream of the TSS.
#' @param both_strands When \code{custom_seq} is specified,
#'     should both strands be considered? TRUE by default.
#' @param anchor The position within the protospacer as determined by
#'     \linkS4class{CrisprNuclease} to use when annotating with overlapping
#'     gene regions.
#' @param n0_max Number of maximum on-target alignments tolerated for
#'     \code{\link{addSpacerAlignmentsIterative}}.
#' @param n1_max Number of maximum 1mm off-target alignments tolerated
#'     for \code{\link{addSpacerAlignmentsIterative}}.
#' @param n2_max Number of maximum 2mm off-target alignments tolerated
#'     for \code{\link{addSpacerAlignmentsIterative}}.
#'
#' @return \code{\link{getSpacerAlignments}} returns a \linkS4class{GRanges} 
#'     object storing spacer alignment data, including genomic coordinates, 
#'     spacer and PAM sequences, and position of mismatches relative to
#'     \code{pam_site}.
#' 
#' @return \code{\link{addSpacerAlignments}} is similar to 
#'     \code{\link{getSpacerAlignments}}, with the addition of adding the 
#'     alignment data to a list-column in \code{mcols(guideSet)} specified
#'     by \code{columnName}. 
#' 
#' @return \code{\link{addSpacerAlignmentsIterative}} is similar to
#'     \code{\link{addSpacerAlignments}}, except that it avoids finding 
#'     alignments for spacer sequences that have a large number of on-targets
#'     and/or off-targets to speed up the off-target search. The parameters
#'     \code{n0_max}, \code{n1_max} and \code{n2_max} specify the maximum
#'     number of on-targets (n0) and off-targets
#'     (n1 for 1-mismatch off-targets, and n2 for 2-mismatch off-targets) 
#'     tolerated before the algorithm stops finding additional off-targets
#'     for spacer sequences that exceed those quotas. 
#' 
#' @details 
#' 
#' The columns stored in \code{mcols(guideSet)[["alignments"]]} are:
#' 
#' \itemize{
#' \item \code{spacer} Spacer sequence of the query gRNA.
#' \item \code{protospacer} Protospacer sequence in the target DNA.
#' \item \code{pam} PAM sequence.
#' \item \code{pam_site} PAM site of the found protospacer.
#' \item \code{n_mismatches} Integer value specifying the number
#'     of nucleotide mismatches between the gRNA spacer sequence 
#'     and the protospacer sequence found in the genome or custom sequence.
#' \item \code{canonical} Whether the PAM sequence of the found protospacer
#'     sequence is canonical.
#' \item \code{cute_site} Cut site of the found protospacer.
#' }
#' 
#' The following columns are also stored when a \code{txObject} is provided:
#' 
#' \itemize{
#' \item \code{cds} Character vector specifying gene names of CDS overlapping
#'     the found protospacer sequence.
#' \item \code{fiveUTRs} Character vector specifying gene names of 5'UTRs
#'     overlapping the found protospacer sequence.
#' \item \code{threeUTRs} Character vector specifying gene names of 3'UTRs
#'     overlapping the found protospacer sequence.
#' \item \code{exons} Character vector specifying gene names of exons
#'     overlapping the found protospacer sequence.
#' \item \code{introns} Character vector specifying gene names of introns
#'     overlapping the found protospacer sequence.
#' \item \code{intergenic} Character vector specifying the nearest gene when
#'     the found protospacer sequence is not located in a gene.
#' \item \code{intergenic_distance} Distance in base pairs from the nearest
#'     gene when the found protospacer sequence is not located in a gene.
#' }
#' 
#' The following columns are also stored when a \code{tssObject} is provided:
#' 
#' \itemize{
#' \item \code{promoters} Character vector specifying gene names of promoters,
#'     as defined by \code{tss_window} relative to the gene TSS, overlapping
#'     the found protospacer sequence.
#' }
#' 
#' @examples 
#' 
#' # Creating a bowtie index:
#' library(Rbowtie)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' fasta <- system.file(package="crisprDesign", "fasta/chr12.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,
#'                       outdir=outdir,
#'                       force=TRUE,
#'                       prefix="chr12")
#' bowtieIndex <- file.path(outdir, "chr12")
#' 
#' # Adding spacer alignments with bowtie:
#' data(guideSetExample, package="crisprDesign")
#' data(grListExample, package="crisprDesign")
#' guideSet <- addSpacerAlignments(guideSetExample,
#'                                 aligner="bowtie",
#'                                 aligner_index=bowtieIndex,
#'                                 bsgenome=BSgenome.Hsapiens.UCSC.hg38,
#'                                 n_mismatches=2,
#'                                 txObject=grListExample)
#' 
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @name addSpacerAlignments
NULL



#' @rdname addSpacerAlignments
#' @export
#' @importFrom S4Vectors mcols
addSpacerAlignmentsIterative <- function(guideSet,
                                         aligner=c("bowtie", "bwa", "biostrings"),
                                         columnName="alignments",
                                         addSummary=TRUE,
                                         txObject=NULL,
                                         tssObject=NULL,
                                         custom_seq=NULL,
                                         aligner_index=NULL,
                                         bsgenome=NULL,
                                         n_mismatches=0,
                                         all_alignments=FALSE,
                                         canonical=TRUE,
                                         ignore_pam=FALSE,
                                         standard_chr_only=TRUE,
                                         both_strands=TRUE,
                                         anchor=c("cut_site", "pam_site"),
                                         tss_window=NULL,
                                         n0_max=5,
                                         n1_max=100,
                                         n2_max=100
){
    guideSet    <- .validateGuideSet(guideSet)
    aligner <- match.arg(aligner)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches, aligner)
    maxAlignments <- c("n0_max", "n1_max", "n2_max")
    maxAlignments <- vapply(maxAlignments, function(x){
        .checkSingleInteger(x, get(x), null_ok=FALSE, sign="non-negative")
        get(x)
    }, FUN.VALUE=numeric(1))
    
    .iterateAddSpacerAlignments <- function(guideSet,
                                            n_mismatches
    ){
        addSpacerAlignments(guideSet,
                            aligner=aligner,
                            columnName=columnName,
                            addSummary=addSummary,
                            txObject=txObject,
                            tssObject=tssObject,
                            custom_seq=custom_seq,
                            aligner_index=aligner_index,
                            bsgenome=bsgenome,
                            n_mismatches=n_mismatches,
                            all_alignments=all_alignments,
                            canonical=canonical,
                            ignore_pam=ignore_pam,
                            standard_chr_only=standard_chr_only,
                            both_strands=both_strands,
                            anchor=anchor,
                            tss_window=tss_window)
    }
    guideSet <- .iterateAddSpacerAlignments(guideSet, 0)
    good <- TRUE
    for (i in seq_len(n_mismatches)){
        mismatch_col <- paste0("n", i-1)
        good <- good &
            S4Vectors::mcols(guideSet)[[mismatch_col]] <= maxAlignments[i]
        guideSet[good] <- .iterateAddSpacerAlignments(guideSet[good], i)
    }
    return(guideSet)
}



#' @rdname addSpacerAlignments
#' @export
#' @importFrom S4Vectors split mcols mcols<-
addSpacerAlignments <- function(guideSet,
                                aligner=c("bowtie", "bwa", "biostrings"),
                                columnName="alignments",
                                addSummary=TRUE,
                                txObject=NULL,
                                tssObject=NULL,
                                custom_seq=NULL,
                                aligner_index=NULL,
                                bsgenome=NULL,
                                n_mismatches=0,
                                n_max_alignments=1000,
                                all_alignments=TRUE,
                                canonical=TRUE,
                                ignore_pam=FALSE,
                                standard_chr_only=TRUE,
                                both_strands=TRUE,
                                anchor=c("cut_site", "pam_site"),
                                tss_window=NULL
){
    guideSet  <- .validateGuideSet(guideSet)
    aligner <- match.arg(aligner)
    .checkString("columnName", columnName)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches, aligner)
    anchor <- match.arg(anchor)
    spacers <- spacers(guideSet, as.character=TRUE)
    uniqueSpacers <- unique(spacers)
    
    aln <- getSpacerAlignments(spacers=uniqueSpacers,
                               aligner=aligner,
                               n_mismatches=n_mismatches,
                               custom_seq=custom_seq,
                               aligner_index=aligner_index,
                               bsgenome=bsgenome,
                               n_max_alignments=n_max_alignments,
                               all_alignments=all_alignments,
                               crisprNuclease=crisprNuclease(guideSet),
                               canonical=canonical,
                               ignore_pam=ignore_pam,
                               standard_chr_only=standard_chr_only,
                               both_strands=both_strands)
    # for both, need to check if txObject/tssObject have compatible seqinfo with bsgenome (or custom_seq...)
    aln <- .addGeneAnnotationColumns(aln,
                                     txObject=txObject,
                                     anchor=anchor)
    aln <- .addPromoterAnnotationColumns(aln,
                                         tssObject=tssObject,
                                         tss_window=tss_window,
                                         anchor=anchor)
    guideSet <- .addAlignmentsSummary(guideSet=guideSet,
                                      aln=aln,
                                      addSummary=addSummary,
                                      n_mismatches=n_mismatches,
                                      spacers=spacers)
    aln <- S4Vectors::split(aln,
                            f=factor(S4Vectors::mcols(aln)$spacer,
                                     levels=uniqueSpacers))
    aln <- aln[spacers]
    names(aln) <- names(guideSet)
    S4Vectors::mcols(guideSet)[[columnName]] <- aln
    return(guideSet)
}



#' @rdname addSpacerAlignments
#' @export
getSpacerAlignments <- function(spacers,
                                aligner=c("bowtie", "bwa", "biostrings"),
                                custom_seq=NULL,
                                aligner_index=NULL,
                                bsgenome=NULL,
                                n_mismatches=0,
                                n_max_alignments=1000,
                                all_alignments=TRUE,
                                crisprNuclease=NULL,
                                canonical=TRUE,
                                ignore_pam=FALSE,
                                standard_chr_only=TRUE,
                                both_strands=TRUE
){
    
    if (.isGuideSet(spacers)){
        spacers <- spacers(spacers)
    }
    spacers <- as.character(spacers)
    aligner <- match.arg(aligner)
    # should obtain from spacers if spacers is GuideSet
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches, aligner)
    
    if (!is.null(custom_seq) & aligner != "biostrings"){
        aligner <- "biostrings"
        message("Setting aligner to 'biostrings' since custom_seq is provided.")
    }
    
    if (aligner %in% c("bowtie", "bwa")){
        aln <- .getSpacerAlignments_indexed(spacers=spacers,
                                            aligner=aligner,
                                            aligner_index=aligner_index,
                                            bsgenome=bsgenome,
                                            n_mismatches=n_mismatches,
                                            n_max_alignments=n_max_alignments,
                                            all_alignments=all_alignments,
                                            crisprNuclease=crisprNuclease,
                                            canonical=canonical,
                                            ignore_pam=ignore_pam,
                                            standard_chr_only=standard_chr_only)
    } else {
        aln <- .getSpacerAlignments_biostrings(spacers=spacers,
                                               custom_seq=custom_seq,
                                               n_mismatches=n_mismatches,
                                               crisprNuclease=crisprNuclease,
                                               canonical=canonical,
                                               ignore_pam=ignore_pam,
                                               both_strands=both_strands)
    }
    return(aln)
}



.validateNumberOfMismatches <- function(n_mismatches,
                                        aligner
){
    errorMessage <- "'n_mismatches' must be an integer value"
    if (aligner == "bowtie"){
        errorMessage <- paste("For bowtie alignments",
                              errorMessage,
                              "between 0 and 3, inclusive")
    }
    if (!is.vector(n_mismatches, mode="numeric")){
        stop(errorMessage)
    }
    isInteger <- n_mismatches == round(n_mismatches)
    isSingleValue <- length(n_mismatches) == 1
    isValidValue <- aligner != "bowtie" || n_mismatches %in% 0:3
    if (!isInteger || !isSingleValue || !isValidValue){
        stop(errorMessage)
    }
    return(n_mismatches)
}



#' @importFrom crisprBowtie runCrisprBowtie
#' @importFrom crisprBwa runCrisprBwa
.getSpacerAlignments_indexed <- function(spacers,
                                         aligner,
                                         aligner_index,
                                         bsgenome,
                                         n_mismatches,
                                         n_max_alignments,
                                         all_alignments,
                                         crisprNuclease,
                                         canonical,
                                         ignore_pam,
                                         standard_chr_only
){
    if (aligner == "bwa" && .Platform$OS.type=="windows"){
        stop("BWA aligner is not available for Windows machines. ",
             "Choose bowtie instead.")   
    }
    spacerLength <- unique(nchar(spacers))
    if (length(spacerLength) > 1){
        stop("All spacer sequences must have the same length.")
    }
    .isBSGenome(bsgenome)
    
    results <- switch(
        aligner,
        "bowtie"=crisprBowtie::runCrisprBowtie(spacers=spacers,
                                               bowtie_index=aligner_index,
                                               bsgenome=bsgenome,
                                               n_mismatches=n_mismatches,
                                               n_max_alignments=n_max_alignments,
                                               crisprNuclease=crisprNuclease,
                                               canonical=canonical,
                                               ignore_pam=ignore_pam,
                                               all_alignments=all_alignments,
                                               force_spacer_length=TRUE),
        "bwa"=crisprBwa::runCrisprBwa(spacers=spacers,
                                      bwa_index=aligner_index,
                                      bsgenome=bsgenome,
                                      n_mismatches=n_mismatches,
                                      crisprNuclease=crisprNuclease,
                                      canonical=canonical,
                                      ignore_pam=ignore_pam,
                                      force_spacer_length=TRUE)
    )
    results <- .alignmentOutput2GRanges(alignments=results,
                                        crisprNuclease=crisprNuclease)
    results <- .setAlignmentSeqInfo(alignments=results,
                                    bsgenome=bsgenome,
                                    standard_chr_only=standard_chr_only)
    alignmentParams <- list(n_mismatches=n_mismatches,
                            canonical=canonical,
                            spacer_len=spacerLength)
    if (aligner == "bowtie"){
        alignmentParams[["n_max_alignments"]] <- n_max_alignments
        alignmentParams[["all_alignments"]] <- all_alignments
    }
    results <- .addAlignmentsMetadata(results,
                                      aligner=aligner,
                                      crisprNuclease=crisprNuclease,
                                      alignmentParams=alignmentParams)
    names(results) <- paste0("aln_", seq_along(results), recycle0=TRUE)
    return(results)
}



#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom Biostrings DNAStringSet
#' @importFrom BiocGenerics strand
.alignmentOutput2GRanges <- function(alignments,
                                     crisprNuclease
){
    mcols <- c("spacer", "protospacer", "pam", "pam_site", "n_mismatches",
               "canonical")
    mcols <- alignments[, mcols, drop=FALSE]
    alignments <- GenomicRanges::GRanges(alignments$chr,
                                         IRanges::IRanges(start=alignments$pam_site,
                                                          width=1),
                                         strand=as.character(alignments$strand))
    S4Vectors::mcols(alignments) <- mcols
    seqCols <- c("spacer", "protospacer", "pam")
    for (i in seqCols){
        S4Vectors::mcols(alignments)[[i]] <-
            Biostrings::DNAStringSet(S4Vectors::mcols(alignments)[[i]])
    }
    strandChar <- as.character(BiocGenerics::strand(alignments))
    alignments$cut_site <- getCutSiteFromPamSite(pam_site=alignments$pam_site,
                                                 strand=strandChar,
                                                 crisprNuclease=crisprNuclease)
    return(alignments)
}



#' @importFrom GenomeInfoDb seqnames seqlevels seqlevels<-
#' @importFrom GenomeInfoDb seqinfo seqinfo<- keepStandardChromosomes
.setAlignmentSeqInfo <- function(alignments,
                                 bsgenome,
                                 standard_chr_only
){
    alignments <- alignments[GenomeInfoDb::seqnames(alignments) != "chrEBV"]
    mtChr <- GenomeInfoDb::seqlevels(alignments) == "chrMT"
    GenomeInfoDb::seqlevels(alignments)[mtChr] <- "chrM"
    
    GenomeInfoDb::seqlevels(alignments) <- GenomeInfoDb::seqlevels(bsgenome)
    GenomeInfoDb::seqinfo(alignments) <- GenomeInfoDb::seqinfo(bsgenome)
    if (standard_chr_only){
        alignments <- GenomeInfoDb::keepStandardChromosomes(alignments,
                                                            pruning.mode="coarse")
    }
    return(alignments)
}



#' @importFrom S4Vectors metadata<-
.addAlignmentsMetadata <- function(alignments,
                                   aligner,
                                   crisprNuclease,
                                   alignmentParams
){
    S4Vectors::metadata(alignments)[["aligner"]] <- aligner
    S4Vectors::metadata(alignments)[["crisprNuclease"]] <- crisprNuclease
    for (i in seq_along(alignmentParams)){
        name <- names(alignmentParams)[i]
        value <- alignmentParams[[i]]
        S4Vectors::metadata(alignments)[[name]] <- value
    }
    return(alignments)
}



#' @importFrom BiocGenerics rbind
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom crisprBase motifs
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @importFrom GenomeInfoDb seqinfo<-
.getSpacerAlignments_biostrings <- function(spacers,
                                            custom_seq,
                                            n_mismatches, 
                                            crisprNuclease,
                                            canonical,
                                            ignore_pam,
                                            both_strands
){
    ## change to allow custom_seq of any length, and as DNAStringSet
    custom_seq <- .validateDNACharacterVariable(seq=custom_seq,
                                                argument="custom_seq",
                                                len=1,
                                                nullOk=FALSE,
                                                exactBases=FALSE)
    
    results <- lapply(spacers, function(x){
        .getCustomSeqAlignments(spacer=x,
                                custom_seq=custom_seq,
                                n_mismatches=n_mismatches,
                                crisprNuclease=crisprNuclease,
                                both_strands=both_strands)
    })
    results <- Reduce(BiocGenerics::rbind, results)
    results <- GenomicRanges::GRanges(
        seqnames=rep("custom", nrow(results)),
        ranges=IRanges::IRanges(start=results$pam_site, width=1),
        strand=results$strand,
        spacer=Biostrings::DNAStringSet(results$spacer),
        protospacer=Biostrings::DNAStringSet(results$seq),
        pam=Biostrings::DNAStringSet(results$pam),
        pam_site=results$pam_site)
    resultsPams <- as.character(S4Vectors::mcols(results)$pam)
    if (!ignore_pam){
        pamMotifs <- crisprBase::motifs(crisprNuclease,
                                        primary=canonical,
                                        expand=TRUE,
                                        as.character=TRUE)
        results <- results[resultsPams %in% pamMotifs]
        resultsPams <- resultsPams[resultsPams %in% pamMotifs]
    }
    results$n_mismatches <- vapply(seq_along(results), function(x){
        adist(S4Vectors::mcols(results)$spacer[x],
              S4Vectors::mcols(results)$protospacer[x])
    }, FUN.VALUE=numeric(1))
    canonicalMotifs <- crisprBase::motifs(crisprNuclease,
                                          primary=TRUE,
                                          expand=TRUE,
                                          as.character=TRUE)
    S4Vectors::mcols(results)$canonical <- resultsPams %in% canonicalMotifs
    S4Vectors::mcols(results)$cut_site <- getCutSiteFromPamSite(
        pam_site=results$pam_site,
        strand=as.character(strand(results)),
        crisprNuclease=crisprNuclease)
    
    GenomeInfoDb::seqinfo(results) <- GenomeInfoDb::Seqinfo(
        seqnames="custom",
        seqlengths=nchar(custom_seq),
        isCircular=FALSE,
        genome="custom")
    
    alignmentParams <- list(n_mismatches=n_mismatches,
                            canonical=canonical,
                            both_strands=both_strands,
                            spacer_len <- unique(nchar(results$spacer)),
                            custom_seq=Biostrings::DNAStringSet(custom_seq))
    results <- .addAlignmentsMetadata(results,
                                      aligner="biostrings",
                                      crisprNuclease=crisprNuclease,
                                      alignmentParams=alignmentParams)
    
    names(results) <- paste0("aln_", seq_along(results), recycle0=TRUE)
    return(results)
}



#' @importFrom BiocGenerics rbind
.getCustomSeqAlignments <- function(spacer,
                                    custom_seq,
                                    n_mismatches,
                                    crisprNuclease,
                                    both_strands
){
    hits <- .getCustomSeqPatternHits(spacer=spacer,
                                     custom_seq=custom_seq,
                                     n_mismatches=n_mismatches,
                                     strand="+")
    if (both_strands){
        hits_rev <- .getCustomSeqPatternHits(spacer=.revComp(spacer),
                                             custom_seq=custom_seq,
                                             n_mismatches=n_mismatches,
                                             strand="-")
        hits <- BiocGenerics::rbind(hits, hits_rev)
    }
    hits$pam_site <- .getPamSiteFromSpacerRange(start=hits$start,
                                                end=hits$end,
                                                strand=hits$strand,
                                                crisprNuclease=crisprNuclease)
    hits$pam <- .getPamFromCustomSeq(custom_seq=custom_seq,
                                     pam_site=hits$pam_site,
                                     strand=hits$strand,
                                     crisprNuclease=crisprNuclease)
    return(hits)
}



#' @importFrom Biostrings matchPattern
#' @importFrom BiocGenerics start end as.data.frame
.getCustomSeqPatternHits <- function(spacer,
                                     custom_seq,
                                     n_mismatches,
                                     strand
){
    hits <-  Biostrings::matchPattern(spacer,
                                      custom_seq,
                                      max.mismatch=n_mismatches)
    inRange5Prime <- BiocGenerics::start(hits) > 0
    inRange3Prime <- BiocGenerics::end(hits) <= nchar(custom_seq)
    hits <- hits[inRange5Prime & inRange3Prime]
    hits <- BiocGenerics::as.data.frame(hits)
    hits$strand <- rep(strand, nrow(hits))
    hits$spacer <- rep(spacer, nrow(hits))
    if (strand == "-" && nrow(hits) > 0){
        hits$seq <- .revComp(hits$seq)
        hits$spacer <- .revComp(spacer)
    }
    return(hits)
}



#' @importFrom crisprBase pamSide pamLength spacerGap
.getPamSiteFromSpacerRange <- function(start,
                                       end,
                                       strand,
                                       crisprNuclease
){
    pamSide <- crisprBase::pamSide(crisprNuclease)
    pamLength <- crisprBase::pamLength(crisprNuclease)
    gap <- crisprBase::spacerGap(crisprNuclease)
    pam_site <- rep(0, length(start))
    
    if (pamSide == "3prime"){
        pam_site[strand == "+"] <- end[strand == "+"] + gap + 1
        pam_site[strand == "-"] <- start[strand == "-"] - gap - 1
    } else {
        pam_site[strand == "+"] <- start[strand == "+"] - gap - pamLength - 1
        pam_site[strand == "-"] <- end[strand == "-"] + gap + pamLength + 1
    }
    return(pam_site)
}



#' @importFrom crisprBase pamLength
.getPamFromCustomSeq <- function(custom_seq,
                                 pam_site,
                                 strand,
                                 crisprNuclease
){
    pamLength <- crisprBase::pamLength(crisprNuclease)
    pams <- vapply(seq_along(pam_site), function(x){
        if (strand[x] == "+"){
            start <- pam_site[x]
            seq <- custom_seq
        } else {
            start <- nchar(custom_seq) - pam_site[x] + 1
            seq <- .revComp(custom_seq)
        }
        substr(seq, start, start+pamLength-1)
    }, FUN.VALUE=character(1))
    return(pams)
}



#' @importFrom S4Vectors mcols mcols<- subjectHits
#' @importFrom GenomicRanges distanceToNearest
.addGeneAnnotationColumns <- function(aln,
                                      txObject,
                                      anchor
){
    if (is.null(txObject)){
        return(aln)
    }
    txObject <- .validateGRangesList(txObject)
    # check that gene_symbol and gene_id are in transcripts/cds/promoters
    
    regions <- c("cds", "fiveUTRs", "threeUTRs", "exons", "introns")
    for (i in regions) {
        regionAnnotation <- .addGeneOverlapByRegion(aln=aln,
                                                    geneRegionModel=txObject[[i]],
                                                    anchor=anchor)
        S4Vectors::mcols(aln)[[i]] <- regionAnnotation
    }
    
    intergenicAnnotation <- lapply(seq_along(aln), function(x){
        geneRegionAnnotation <- S4Vectors::mcols(aln)[x, regions]
        if (any(!is.na(geneRegionAnnotation))){
            results <- data.frame(gene=NA_character_,
                                  dist=NA_integer_)
            return(results)
        }
        hit <- GenomicRanges::distanceToNearest(aln[x],
                                                txObject[["transcripts"]],
                                                ignore.strand=TRUE)
        nearestGene <- txObject[["transcripts"]][S4Vectors::subjectHits(hit)]
        geneSymbol <- S4Vectors::mcols(nearestGene)[["gene_symbol"]]
        if (is.na(geneSymbol) || geneSymbol == ""){
            nearestGene <- S4Vectors::mcols(nearestGene)[["gene_id"]]
        } else {
            nearestGene <- geneSymbol
        }
        nearestGene <- paste0(nearestGene, collapse=";")
        nearestDistance <- S4Vectors::mcols(hit)[["distance"]]
        results <- data.frame(gene=nearestGene,
                              dist=nearestDistance)
        return(results)
    })
    intergenicAnnotation <- Reduce(rbind, intergenicAnnotation)
    S4Vectors::mcols(aln)[["intergenic"]] <- intergenicAnnotation$gene
    S4Vectors::mcols(aln)[["intergenic_distance"]] <- intergenicAnnotation$dist
    
    return(aln)
}



#' @importFrom IRanges IRanges ranges<-
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols
.addGeneOverlapByRegion <- function(aln,
                                    geneRegionModel,
                                    anchor
){
    anchor <- .validateAnchor(anchor, aln)
    IRanges::ranges(aln) <- IRanges::IRanges(start=S4Vectors::mcols(aln)[[anchor]],
                                             width=1)
    overlaps <- suppressWarnings(
        GenomicRanges::findOverlaps(aln,
                                    geneRegionModel,
                                    ignore.strand=TRUE)
        )
    
    regionAnnotation <- vapply(seq_along(aln), function(x){
        if (!x %in% S4Vectors::queryHits(overlaps)){
            return(NA_character_)
        }
        indicesOfHits <- S4Vectors::queryHits(overlaps) == x
        regionHits <- S4Vectors::subjectHits(overlaps)[indicesOfHits]
        geneRegionModelSubset <- S4Vectors::mcols(geneRegionModel)[regionHits,]
        geneHits <- geneRegionModelSubset[["gene_symbol"]]      ## can be missing
        geneIds <- geneRegionModelSubset[["gene_id"]]           ## required
        missingSymbols <- is.na(geneHits) | geneHits == ""
        geneHits[missingSymbols] <- geneIds[missingSymbols]
        if (all(geneHits == "")){
            return(NA_character_)
        }
        geneHits <- geneHits[geneHits != ""]
        geneHits <- unique(geneHits)
        paste0(geneHits, collapse=";")
    }, FUN.VALUE=character(1))
    
    return(regionAnnotation)
}



#' @importFrom GenomicRanges promoters
#' @importFrom S4Vectors mcols<-
.addPromoterAnnotationColumns <- function(aln,
                                          tssObject,
                                          tss_window,
                                          anchor
){
    if (is.null(tssObject)){
        return(aln)
    }
    tssObject <- .validateTssObject(tssObject)
    tss_window <- .validateTssWindow(tss_window)
    
    tssObject <- GenomicRanges::promoters(tssObject,
                                          upstream=(-1*tss_window[1]),
                                          downstream=tss_window[2])
    promoterAnnotation <- .addGeneOverlapByRegion(aln,
                                                  geneRegionModel=tssObject,
                                                  anchor=anchor)
    S4Vectors::mcols(aln)[["promoters"]] <- promoterAnnotation
    return(aln)
}



#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics cbind
.addAlignmentsSummary <- function(guideSet,
                                  aln,
                                  addSummary,
                                  spacers,
                                  n_mismatches
){
    if (addSummary){
        hasAlnSummaryCols <- any(grepl("^n[0-9](_[cp])?$",
                                       colnames(S4Vectors::mcols(guideSet))))
        if (hasAlnSummaryCols){
            warning("Overwriting existing alignments summary. ",
                    "To avoid overwriting, set addSummary=FALSE.")
        }
        
        alignmentSummary <- .getAlignmentsSummary(aln=aln,
                                                  spacers=spacers,
                                                  n_mismatches=n_mismatches)
        S4Vectors::mcols(guideSet) <- BiocGenerics::cbind(
            S4Vectors::mcols(guideSet),
            alignmentSummary)
    }
    return(guideSet)
}



#' @importFrom S4Vectors mcols
.getAlignmentsSummary <- function(aln,
                                  spacers,
                                  n_mismatches
){
    seq_mismatches <- c(0, seq_len(n_mismatches))
    alignmentsSummary <- .tallyAlignments(aln=aln,
                                          seq_mismatches=seq_mismatches,
                                          spacers=spacers,
                                          suffix="")
    
    hasGeneAnnotation <- "cds" %in% colnames(S4Vectors::mcols(aln))
    if (hasGeneAnnotation){
        targetInCds <- !is.na(S4Vectors::mcols(aln)[["cds"]])
        alnCds <- aln[targetInCds]
        cdsAlignments <- .tallyAlignments(aln=alnCds,
                                          seq_mismatches=seq_mismatches,
                                          spacers=spacers,
                                          suffix="_c")
        alignmentsSummary <- cbind(alignmentsSummary, cdsAlignments)
    }
    
    hasPromoterAnnotation <- "promoters" %in% colnames(S4Vectors::mcols(aln))
    if (hasPromoterAnnotation){
        targetInPromoter <- !is.na(S4Vectors::mcols(aln)[["promoters"]])
        alnPromoter <- aln[targetInPromoter]
        promoterAlignments <- .tallyAlignments(aln=alnPromoter,
                                               seq_mismatches=seq_mismatches,
                                               spacers=spacers,
                                               suffix="_p")
        alignmentsSummary <- cbind(alignmentsSummary, promoterAlignments)
    }
    
    return(alignmentsSummary)
}



#' @importFrom S4Vectors mcols
.tallyAlignments <- function(aln,
                             seq_mismatches,
                             spacers,
                             suffix
){
    tally <- lapply(seq_mismatches, function(i){
        tallyByMismatchType <- vapply(spacers, function(ii){
            mismatchHits <- S4Vectors::mcols(aln)[["n_mismatches"]] == i
            spacerHits <- S4Vectors::mcols(aln)[["spacer"]] == ii
            sum(mismatchHits & spacerHits)
        }, FUN.VALUE=numeric(1))
        data.frame(tallyByMismatchType)
    })
    tally <- Reduce(cbind, tally)
    colnames(tally) <- paste0("n", seq_mismatches, suffix)
    return(tally)
}
