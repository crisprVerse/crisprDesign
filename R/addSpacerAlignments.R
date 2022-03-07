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
#' @param colname String specifying the columm name storing the alignments
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
#' @param canonical \code{TRUE} returns only those alignments having canonical
#'     PAM sequences; \code{FALSE} returns alignments having canonical or
#'     noncanonical PAM sequences; \code{NA} returns all alignments regardless
#'     of their PAM sequence.
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
#' @param annotationType Gene identifier to return when annotating alignments
#'     with gene and/or promoter overlaps. Corresponding \code{txObject} or
#'     \code{tssObject} argument must have mcol column name for selected type.
#' @param alignmentThresholds Named numeric vector of the maximum on-target
#'     alignments tolerated for \code{\link{addSpacerAlignmentsIterative}}.
#'     Thresholds not provided will take default values.
#'
#' @return \code{\link{getSpacerAlignments}} returns a \linkS4class{GRanges} 
#'     object storing spacer alignment data, including genomic coordinates, 
#'     spacer and PAM sequences, and position of mismatches relative to
#'     \code{pam_site}.
#' 
#' @return \code{\link{addSpacerAlignments}} is similar to 
#'     \code{\link{getSpacerAlignments}}, with the addition of adding the 
#'     alignment data to a list-column in \code{mcols(guideSet)} specified
#'     by \code{colname}. 
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
#' \dontrun{
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
#' }
#' 
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @name addSpacerAlignments
NULL



#' @rdname addSpacerAlignments
#' @export
#' @importFrom S4Vectors mcols mcols<-
addSpacerAlignmentsIterative <- function(guideSet,
                                         aligner=c("bowtie", "bwa", "biostrings"),
                                         colname="alignments",
                                         addSummary=TRUE,
                                         txObject=NULL,
                                         tssObject=NULL,
                                         custom_seq=NULL,
                                         aligner_index=NULL,
                                         bsgenome=NULL,
                                         n_mismatches=0,
                                         all_alignments=FALSE,
                                         canonical=TRUE,
                                         standard_chr_only=TRUE,
                                         both_strands=TRUE,
                                         anchor=c("cut_site", "pam_site"),
                                         annotationType=c("gene_symbol", "gene_id"),
                                         tss_window=NULL,
                                         alignmentThresholds=c(n0=5,
                                                               n1=100,
                                                               n2=100,
                                                               n3=1000,
                                                               n4=1000)
){
    aligner <- match.arg(aligner)
    if (aligner=="bowtie" & n_mismatches>3){
        stop("For bowtie aligner, n_mismatches must be either 0,1,2 or 3.")
    }
    guideSet    <- .validateGuideSet(guideSet)
    aligner <- match.arg(aligner)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches, aligner)
    maxAlignments <- .validateAlignmentThresholds(alignmentThresholds)
    
    .iterateAddSpacerAlignments <- function(guideSet,
                                            n_mismatches
    ){
        suppressWarnings(
            addSpacerAlignments(guideSet,
                                aligner=aligner,
                                colname=colname,
                                addSummary=addSummary,
                                txObject=txObject,
                                tssObject=tssObject,
                                custom_seq=custom_seq,
                                aligner_index=aligner_index,
                                bsgenome=bsgenome,
                                n_mismatches=n_mismatches,
                                all_alignments=all_alignments,
                                canonical=canonical,
                                standard_chr_only=standard_chr_only,
                                both_strands=both_strands,
                                anchor=anchor,
                                annotationType=annotationType,
                                tss_window=tss_window)
        )
    }
    guideSet <- .iterateAddSpacerAlignments(guideSet, 0)
    good <- TRUE
    for (i in seq_len(n_mismatches)){
        mismatch_col <- paste0("n", i-1)
        good <- good &
            S4Vectors::mcols(guideSet)[[mismatch_col]] <= maxAlignments[i]
        if (any(good)){
            updatedGuideSet <- .iterateAddSpacerAlignments(guideSet[good], i)
            newCols <- setdiff(colnames(S4Vectors::mcols(updatedGuideSet)),
                               colnames(S4Vectors::mcols(guideSet)))
            for (ii in newCols){
                S4Vectors::mcols(guideSet)[[ii]] <- as.numeric(NA)
            }
            guideSet[good] <- updatedGuideSet
        } else {
            na_col <- paste0("n", i)
            S4Vectors::mcols(guideSet)[[na_col]] <- as.numeric(NA)
        }
    }
    aln <- S4Vectors::mcols(guideSet)[[colname]]
    S4Vectors::mcols(guideSet)[[colname]] <- NULL
    S4Vectors::mcols(guideSet)[[colname]] <- aln
    return(guideSet)
}




.validateAlignmentThresholds <- function(alignmentThresholds
){
    maxAlignments <- c(n0=5, n1=100, n2=100, n3=1000, n4=1000)
    if (is.null(alignmentThresholds)){
        return(maxAlignments)
    }
    if (is.null(names(alignmentThresholds)) ||
        !is.vector(alignmentThresholds, mode="numeric")){
        stop("alignmentThresholds must be a named numeric vector")
    }
    duplicatedNames <- duplicated(names(alignmentThresholds))
    invalidNames <- !names(alignmentThresholds) %in% names(maxAlignments)
    if (any(duplicatedNames) || any(invalidNames)){
        stop(paste("names for alignmentThresholds must be unique and in",
                   "c('n0', 'n1', 'n2', 'n3', 'n4')"))
    }
    
    maxAlignments <- maxAlignments[setdiff(names(maxAlignments),
                                           names(alignmentThresholds))]
    alignmentThresholds <- c(alignmentThresholds, maxAlignments)
    alignmentThresholds <- alignmentThresholds[order(names(alignmentThresholds))]
    
    isInteger <- all(alignmentThresholds == round(alignmentThresholds))
    isNonNegative <- all(alignmentThresholds >= 0)
    isNotNA <- all(!is.na(alignmentThresholds))
    if (!isInteger | !isNonNegative | !isNotNA){
        stop("alignmentThresholds must be a named numeric vector of non-negative integers")
    }

    return(alignmentThresholds)
}



#' @rdname addSpacerAlignments
#' @export
#' @importFrom S4Vectors split mcols mcols<-
addSpacerAlignments <- function(guideSet,
                                aligner=c("bowtie", "bwa", "biostrings"),
                                colname="alignments",
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
                                standard_chr_only=TRUE,
                                both_strands=TRUE,
                                anchor=c("cut_site", "pam_site"),
                                annotationType=c("gene_symbol", "gene_id"),
                                tss_window=NULL
){
    guideSet  <- .validateGuideSet(guideSet)
    aligner <- match.arg(aligner)
    .checkString("colname", colname)
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
                               standard_chr_only=standard_chr_only,
                               both_strands=both_strands)
    crisprNuclease <- crisprNuclease(guideSet)
    if (aligner != "biostrings" & isDnase(crisprNuclease)){
        annotationType <- match.arg(annotationType)
        aln <- .addGeneAnnotationColumns(aln,
                                         txObject=txObject,
                                         anchor=anchor,
                                         annotationType=annotationType)
        aln <- .addPromoterAnnotationColumns(aln,
                                             tssObject=tssObject,
                                             tss_window=tss_window,
                                             anchor=anchor,
                                             annotationType=annotationType)
    }
    if (isRnase(crisprNuclease)){
        aln <- .addTranscriptAnnotationColumns(aln,
                                               txObject=txObject)
    }
    
         
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
    S4Vectors::mcols(guideSet)[[colname]] <- aln
    return(guideSet)
}



#' @rdname addSpacerAlignments
#' @export
#' @importFrom methods is
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
                                standard_chr_only=TRUE,
                                both_strands=TRUE
){
    
    if (.isGuideSet(spacers)){
        spacers <- spacers(spacers)
    }
    if (methods::is(spacers, "XStringSet") || methods::is(spacers, "XString")){
        spacers <- as.character(spacers)
    }
    if (!is.vector(spacers, mode="character")){
        stop(paste("'spacers' argument must be a GuideSet, XString, or",
                   "XStringSet object, or a character vector"))
    }
    aligner <- match.arg(aligner)
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches, aligner)
    stopifnot("canonical must be either TRUE, FALSE, or NA" = {
        is.logical(canonical) && length(canonical) == 1
    })
    
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
                                            standard_chr_only=standard_chr_only)
    } else {
        aln <- .getSpacerAlignments_biostrings(spacers=spacers,
                                               custom_seq=custom_seq,
                                               n_mismatches=n_mismatches,
                                               crisprNuclease=crisprNuclease,
                                               canonical=canonical,
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
    if (!isRnase(crisprNuclease)){
        .isBSGenome(bsgenome)
    }
    
    
    if (isRnase(crisprNuclease)){
        bowtie_mode <- "protospacer"
    } else {
        bowtie_mode <- "spacer"
    }

    results <- switch(
        aligner,
        "bowtie"=crisprBowtie::runCrisprBowtie(spacers=spacers,
                                               bowtie_index=aligner_index,
                                               mode=bowtie_mode,
                                               bsgenome=bsgenome,
                                               n_mismatches=n_mismatches,
                                               n_max_alignments=n_max_alignments,
                                               crisprNuclease=crisprNuclease,
                                               canonical=canonical,
                                               ignore_pam=is.na(canonical),
                                               all_alignments=all_alignments,
                                               force_spacer_length=TRUE),
        "bwa"=crisprBwa::runCrisprBwa(spacers=spacers,
                                      bwa_index=aligner_index,
                                      bsgenome=bsgenome,
                                      n_mismatches=n_mismatches,
                                      crisprNuclease=crisprNuclease,
                                      canonical=canonical,
                                      ignore_pam=is.na(canonical),
                                      force_spacer_length=TRUE)
    )
    results <- .alignmentOutput2GRanges(alignments=results,
                                        crisprNuclease=crisprNuclease)
    if (!isRnase(crisprNuclease)){
        results <- .setAlignmentSeqInfo(alignments=results,
                                        bsgenome=bsgenome,
                                        standard_chr_only=standard_chr_only)
    }
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
    .checkBoolean("standard_chr_only", standard_chr_only)
    
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
                                            both_strands
){
    custom_seq <- .setCustomSeqNames(custom_seq)
    .checkBoolean("both_strands", both_strands)
    results <- lapply(spacers, function(x){
        .getCustomSeqAlignments(spacer=x,
                                custom_seq=custom_seq,
                                n_mismatches=n_mismatches,
                                crisprNuclease=crisprNuclease,
                                both_strands=both_strands)
    })
    results <- Reduce(BiocGenerics::rbind, results)
    results <- GenomicRanges::GRanges(
        seqnames=results$seqnames,
        ranges=IRanges::IRanges(start=results$pam_site, width=1), # handle null case
        strand=results$strand,
        spacer=Biostrings::DNAStringSet(results$spacer),
        protospacer=Biostrings::DNAStringSet(results$seq),
        pam=Biostrings::DNAStringSet(results$pam),
        pam_site=results$pam_site)
    resultsPams <- as.character(S4Vectors::mcols(results)$pam)
    if (!is.na(canonical)){
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
    
    GenomeInfoDb::seqlevels(results) <- names(custom_seq)
    GenomeInfoDb::seqinfo(results) <- GenomeInfoDb::Seqinfo(
        seqnames=names(custom_seq),
        seqlengths=nchar(custom_seq),
        isCircular=rep(FALSE, length(custom_seq)),
        genome="custom")
    
    alignmentParams <- list(n_mismatches=n_mismatches,
                            canonical=canonical,
                            both_strands=both_strands,
                            spacer_len=unique(nchar(results$spacer)),
                            custom_seq=Biostrings::DNAStringSet(custom_seq))
    results <- .addAlignmentsMetadata(results,
                                      aligner="biostrings",
                                      crisprNuclease=crisprNuclease,
                                      alignmentParams=alignmentParams)
    
    names(results) <- paste0("aln_", seq_along(results), recycle0=TRUE)
    return(results)
}



.setCustomSeqNames <- function(custom_seq
){
    custom_seq_names <- names(custom_seq)
    custom_seq <- .validateDNACharacterVariable(seq=custom_seq,
                                                argument="custom_seq",
                                                len=NULL,
                                                nullOk=FALSE,
                                                exactBases=FALSE)
    if (is.null(custom_seq_names)){
        names(custom_seq) <- paste0("custom_seq", seq_along(custom_seq))
    }
    missingNames <- which(is.na(custom_seq_names) | custom_seq_names=="")
    newNames <- paste0("custom_seq", missingNames, recycle0=TRUE)
    names(custom_seq)[missingNames] <- newNames
    return(custom_seq)
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
    hits <- .addPamSiteFromSpacerRange(hits=hits,
                                       crisprNuclease=crisprNuclease)
    hits <- .addPamFromCustomSeq(hits=hits,
                                 custom_seq=custom_seq,
                                 crisprNuclease=crisprNuclease)
    return(hits)
}



#' @importFrom Biostrings matchPattern
#' @importFrom BiocGenerics as.data.frame
.getCustomSeqPatternHits <- function(spacer,
                                     custom_seq,
                                     n_mismatches,
                                     strand
){
    hits <-  Biostrings::vmatchPattern(spacer,
                                       custom_seq,
                                       max.mismatch=n_mismatches)
    hits <- lapply(seq_along(hits), function(x){
        inRange5Prime <- BiocGenerics::start(hits[[x]]) > 0
        inRange3Prime <- BiocGenerics::end(hits[[x]]) <= nchar(custom_seq[x])
        perSeqHits <- hits[[x]][inRange5Prime & inRange3Prime]
        perSeqHits <- BiocGenerics::as.data.frame(perSeqHits)
        perSeqHits$seqnames <- rep(names(custom_seq)[x], nrow(perSeqHits))
        perSeqHits
    })
    hits <- Reduce(rbind, hits)
    # seq is protospacer in custom_seq
    hits$seq <- vapply(seq_len(nrow(hits)), function(x){
        sourceSeq <- hits$seqnames[x]
        sourceSeq <- custom_seq[[sourceSeq]]
        substr(sourceSeq, hits$start[x], hits$end[x])
    }, FUN.VALUE=character(1))
    hits$strand <- rep(strand, nrow(hits))
    hits$spacer <- rep(spacer, nrow(hits))
    if (strand == "-" && nrow(hits) > 0){
        hits$seq <- .revComp(hits$seq)
        hits$spacer <- .revComp(spacer)
    }
    return(hits)
}




#' @importFrom crisprBase pamSide pamLength spacerGap
.addPamSiteFromSpacerRange <- function(hits,
                                       crisprNuclease
){
    pamSide <- crisprBase::pamSide(crisprNuclease)
    pamLength <- crisprBase::pamLength(crisprNuclease)
    gap <- crisprBase::spacerGap(crisprNuclease)
    pam_site <- rep(0, nrow(hits))
    forwardStrand <- hits$strand == "+"
    reverseStrand <- hits$strand == "-"
    
    if (pamSide == "3prime"){
        pam_site[forwardStrand] <- hits$end[forwardStrand] + gap + 1
        pam_site[reverseStrand] <- hits$start[reverseStrand] - gap - 1
    } else {
        pam_site[forwardStrand] <- hits$start[forwardStrand] - gap - pamLength - 1
        pam_site[reverseStrand] <- hits$end[reverseStrand] + gap + pamLength + 1
    }
    hits$pam_site <- pam_site
    return(hits)
}




#' @importFrom crisprBase pamLength
.addPamFromCustomSeq <- function(hits,
                                 custom_seq,
                                 crisprNuclease
){
    pamLength <- crisprBase::pamLength(crisprNuclease)
    pams <- vapply(seq_len(nrow(hits)), function(x){
        seq <- custom_seq[[hits$seqnames[x]]]
        start <- hits$pam_site[x]
        if (hits$strand[x] == "+"){
            end <- start + pamLength - 1
            substr(seq, start, end)
        } else {
            end <- start - pamLength + 1
            .revComp(substr(seq, end, start))
        }
    }, FUN.VALUE=character(1))
    hits$pam <- pams
    return(hits)
}



# Function to add transcript annotation
# for RNases
.addTranscriptAnnotationColumns <- function(aln,
                                            txObject=txObject
){  
    if (is.null(txObject)){
        return(aln)
    }
    tx2GeneTable <- .getTx2GeneTable(txObject)
    txids <- as.character(seqnames(aln))
    if (!any(txids %in% tx2GeneTable$tx_id)){
        warning("None of the transcripts in the alignment",
                " table is found in the txObject.")
    }
    wh <- match(txids, tx2GeneTable$tx_id)
    aln$gene_id <- tx2GeneTable[wh, "gene_id"]
    aln$gene_symbol <- tx2GeneTable[wh, "gene_symbol"]
    return(aln)
}





#' @importFrom GenomeInfoDb checkCompatibleSeqinfo
#' @importFrom S4Vectors mcols<-
.addGeneAnnotationColumns <- function(aln,
                                      txObject,
                                      anchor,
                                      annotationType
){
    if (is.null(txObject)){
        return(aln)
    }
    txObject <- .validateGRangesList(txObject)
    GenomeInfoDb::checkCompatibleSeqinfo(aln, txObject)
    
    regions <- c("cds", "fiveUTRs", "threeUTRs", "exons", "introns")
    for (i in regions) {
        regionAnnotation <- .addGeneOverlapByRegion(aln=aln,
                                                    geneRegionModel=txObject[[i]],
                                                    anchor=anchor,
                                                    annotationType=annotationType)
        S4Vectors::mcols(aln)[[i]] <- regionAnnotation
    }
    aln <- .addIntergenicAnnotation(aln=aln,
                                    txModel=txObject[["transcripts"]],
                                    anchor=anchor,
                                    annotationType=annotationType)
    return(aln)
}



#' @importFrom IRanges IRanges ranges<-
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols
.addGeneOverlapByRegion <- function(aln,
                                    geneRegionModel,
                                    anchor,
                                    annotationType
){
    anchor <- .validateAnchor(anchor, aln)
    IRanges::ranges(aln) <- IRanges::IRanges(start=S4Vectors::mcols(aln)[[anchor]],
                                             width=1)
    if (!annotationType %in% colnames(S4Vectors::mcols(geneRegionModel))){
        stop(sprintf("'%s' not found in gene model", annotationType))
    }
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
        geneHits <- geneRegionModelSubset[[annotationType]]
        geneHits <- geneHits[geneHits != ""]
        geneHits <- unique(geneHits)
        paste0(geneHits, collapse=";")
    }, FUN.VALUE=character(1))
    
    return(regionAnnotation)
}



#' @importFrom S4Vectors mcols mcols<- subjectHits queryHits
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges distanceToNearest
.addIntergenicAnnotation <- function(aln,
                                     txModel,
                                     anchor,
                                     annotationType
){
    anchor <- .validateAnchor(anchor, aln)
    anchorSite <- S4Vectors::mcols(aln)[[anchor]]
    alnSetToAnchor <- aln
    IRanges::ranges(alnSetToAnchor) <- IRanges::IRanges(start=anchorSite,
                                                        width=1)
    nearestGene <- GenomicRanges::distanceToNearest(alnSetToAnchor,
                                                    txModel,
                                                    ignore.strand=TRUE)
    nearestGene <- filterOutAlnWithGeneRegionAnnotation(aln, nearestGene)
    
    intergenic <- S4Vectors::mcols(txModel)[[annotationType]]
    intergenic <- intergenic[S4Vectors::subjectHits(nearestGene)]
    intergenic_distance <- S4Vectors::mcols(nearestGene)[["distance"]]
    nearestGene <- data.frame(aln_index=S4Vectors::queryHits(nearestGene),
                              intergenic=intergenic,
                              intergenic_distance=intergenic_distance)

    S4Vectors::mcols(aln)[["intergenic"]] <- NA_character_
    S4Vectors::mcols(aln)[["intergenic_distance"]] <- NA_integer_
    for (i in c("intergenic", "intergenic_distance")){
        S4Vectors::mcols(aln)[[i]][nearestGene$aln_index] <- nearestGene[[i]]
    }
    
    return(aln)
}


#' @importFrom S4Vectors mcols
filterOutAlnWithGeneRegionAnnotation <- function(aln,
                                                 nearestGene
){
    geneRegions <- c("cds", "fiveUTRs", "threeUTRs", "exons", "introns")
    geneCoverage <- S4Vectors::mcols(aln)[geneRegions]
    hasGeneAnnotation <- apply(geneCoverage, 1, function(x){
        all(is.na(x))
    })
    nearestGene <- nearestGene[hasGeneAnnotation]
    return(nearestGene)
}



#' @importFrom GenomeInfoDb checkCompatibleSeqinfo
#' @importFrom GenomicRanges promoters
#' @importFrom S4Vectors mcols<-
.addPromoterAnnotationColumns <- function(aln,
                                          tssObject,
                                          tss_window,
                                          anchor,
                                          annotationType
){
    if (is.null(tssObject)){
        return(aln)
    }
    tssObject <- .validateTssObject(tssObject)
    GenomeInfoDb::checkCompatibleSeqinfo(aln, tssObject)
    tss_window <- .validateTssWindow(tss_window)
    
    tssObject <- GenomicRanges::promoters(tssObject,
                                          upstream=(-1*tss_window[1]),
                                          downstream=tss_window[2])
    promoterAnnotation <- .addGeneOverlapByRegion(aln,
                                                  geneRegionModel=tssObject,
                                                  anchor=anchor,
                                                  annotationType=annotationType)
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
        hasAlnSummaryCols <- any(grepl("^n[0-9](_[cptxgen])?$",
                                       colnames(S4Vectors::mcols(guideSet))))
        if (hasAlnSummaryCols){
            warning("Overwriting existing alignments summary. ",
                    "To avoid overwriting, set addSummary=FALSE.")
        }
        
        if (isDnase(crisprNuclease(guideSet))){
            alignmentSummary <- .getAlignmentsSummary_dnase(aln=aln,
                                                            spacers=spacers,
                                                            n_mismatches=n_mismatches)
        } else {
            alignmentSummary <- .getAlignmentsSummary_rnase(aln=aln,
                                                            spacers=spacers,
                                                            n_mismatches=n_mismatches)
        }
       
        for (i in colnames(alignmentSummary)){
            S4Vectors::mcols(guideSet)[[i]] <- alignmentSummary[[i]]
        }
        # S4Vectors::mcols(guideSet) <- BiocGenerics::cbind(
        #     S4Vectors::mcols(guideSet),
        #     alignmentSummary)
    }
    return(guideSet)
}



#' @importFrom S4Vectors mcols
.getAlignmentsSummary_dnase <- function(aln,
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


.getAlignmentsSummary_rnase <- function(aln,
                                        spacers,
                                        n_mismatches
){
    seq_mismatches <- c(0, seq_len(n_mismatches))
    alignmentsSummary <- .tallyAlignments(aln=aln,
                                          seq_mismatches=seq_mismatches,
                                          spacers=spacers,
                                          suffix="_tx")
    hasGeneAnnotation <- "gene_id" %in% colnames(S4Vectors::mcols(aln))
    if (hasGeneAnnotation){
        geneAlignments <- .tallyAlignments(aln=aln,
                                           seq_mismatches=seq_mismatches,
                                           spacers=spacers,
                                           groupBy="gene_id",
                                           suffix="_gene")
        alignmentsSummary <- cbind(alignmentsSummary, geneAlignments)
    }
    return(alignmentsSummary)
}


#' @importFrom S4Vectors mcols
.tallyAlignments <- function(aln,
                             seq_mismatches,
                             spacers,
                             groupBy=NULL,
                             suffix
){
    tally <- lapply(seq_mismatches, function(i){
        tallyByMismatchType <- vapply(spacers, function(ii){
            mismatchHits <- S4Vectors::mcols(aln)[["n_mismatches"]] == i
            spacerHits <- S4Vectors::mcols(aln)[["spacer"]] == ii
            if (is.null(groupBy)){
                out <- sum(mismatchHits & spacerHits)
            } else {
                genes <- S4Vectors::mcols(aln)[["gene_id"]]
                out <- length(unique(genes[mismatchHits & spacerHits]))
            }
            return(out)
        }, FUN.VALUE=numeric(1))
        data.frame(tallyByMismatchType)
    })
    tally <- Reduce(cbind, tally)
    colnames(tally) <- paste0("n", seq_mismatches, suffix)
    return(tally)
}
