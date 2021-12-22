#' @title Functions for finding and characterizing on- and off-targets of
#'     spacer sequences.
#' 
#' @description Functions for finding and characterizing on- and off-targets of
#'     spacer sequences.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param aligner Which genomic alignment method should be used?
#'     Must be one of "bowtie", "bwa", and "biostrings".
#'    "bowtie" by default.
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
#' @param bowtie_index String specifying bowtie index.
#'     Must be provided when \code{aligner=="bowtie"}.
#' @param bwa_index String specifying BWA index.
#'     Must be provided when \code{aligner=="bwa"}.
#' @param seqlevelsStyle String specifying which type of seqnames
#'     should be used. Default is "UCSC" (e.g. "chr7"; "NCBI"
#'     style would be "7").
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
#' @param cut_offset Distance in nucleotides between \code{pam_site}
#'     and \code{cut_site}, if different from default offset specified in
#'     the \code{crisprNuclease} object. 
#' @param standard_chr_only Should only standard chromosomes be considered?
#'     TRUE by default.
#' @param tss_window Window size of promoters upstream of gene TSS to search
#'     for overlap with spacer sequence. Must be a numeric vector of length 2:
#'     upstream limit and downstream limit. Default is \code{c(-500, 500)},
#'     which includes 500bp upstream and downstream of the TSS.
#' @param both_strands When \code{custom_seq} is specified,
#'     should both strands be considered? TRUE by default.
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
#' The differential columns stored in
#' \code{mcols(guideSet)[["geneAnnotation"]]} are:
#' 
#' \itemize{
#' \item \code{spacer} Transcript ID.
#' \item \code{pam}
#' \item \code{pam_site} PAM site of the found protospacer.
#' \item \code{cute_site} Cut site of the found protospacer.
#' \item \code{n_mismatches} Integer value specifying the number
#'     of nucleotide mismatches between the spacer sequence of the query
#'     and the spacer sequence found in the genome or custom sequence.
#' \item \code{mm1} Numeric value specifying the relative position of the
#'     first mismatch with respect to the PAM site. NA if there is no
#'     mismatch.
#' \item \code{mm2} Numeric value specifying the relative position of the
#'     second mismatch with respect to the PAM site. NA if there is no
#'     second mismatch.
#' \item \code{mm3} Numeric value specifying the relative position of the
#'     third mismatch with respect to the PAM site. NA if there is no
#'     third mismatch.
#' \item \code{canonical} Is the PAM sequence of the found protospacer sequence
#'     canonical?
#' \item \code{query} Spacer sequence of the query gRNA.
#' \item \code{cds} Character vector specifying gene names of CDS overlapping
#'     the found protospacer sequence.
#' \item \code{fiveUTRs} Character vector specifying gene names of 5'UTRs overlapping
#'     the found protospacer sequence.
#' \item \code{threeUTRs} Character vector specifying gene names of 3'UTRs overlapping
#'     the found protospacer sequence.
#' \item \code{exons} Character vector specifying gene names of exons overlapping
#'     the found protospacer sequence.
#' \item \code{introns} Character vector specifying gene names of introns overlapping
#'     the found protospacer sequence.
#' \item \code{intergenic} Charater vector specifying neighborhood genes when the found
#'     protospacer sequence is not located in a gene.
#' }
#' 
#' @examples 
#' 
#' 
#' # Creating a bowtie index:
#' library(Rbowtie)
#' fasta <- system.file(package="crisprDesign", "fasta/chr12.fa")
#' outdir <- tempdir()
#' Rbowtie::bowtie_build(fasta,
#'                       outdir=outdir,
#'                       force=TRUE,
#'                       prefix="chr12")
#' bowtieIndex <- file.path(outdir, "chr12")
#' 
#' # Adding spacer alignments with bowtie:
#' guideSet <- addSpacerAlignments(guideSetExample,
#'                                 aligner="bowtie",
#'                                 bowtie_index=bowtieIndex,
#'                                 n_mismatches=2,
#'                                 txObject=grListExample)
#' 
#' # Creating a bwa index:
#' library(Rbwa)
#' bwaIndex <- file.path(outdir, "chr12")
#' Rbwa::bwa_build_index(fasta,
#'                       index_prefix=bwaIndex)
#' 
#' # Adding spacer alignments with bowtie:
#' guideSet <- addSpacerAlignments(guideSetExample,
#'                                 aligner="bwa",
#'                                 bwa_index=bwaIndex,
#'                                 n_mismatches=2,
#'                                 txObject=grListExample)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @name addSpacerAlignments
NULL




#' @rdname addSpacerAlignments
#' @export
addSpacerAlignmentsIterative <- function(guideSet,
                                         aligner=c("bowtie", "bwa", "biostrings"),
                                         columnName="alignments",
                                         addSummary=TRUE,
                                         txObject=NULL,
                                         tssObject=NULL,
                                         custom_seq=NULL,
                                         bowtie_index=NULL,
                                         bwa_index=NULL,
                                         seqlevelsStyle=c("UCSC", "NCBI"),
                                         bsgenome=NULL,
                                         n_mismatches=0,
                                         all_alignments=FALSE,
                                         canonical=TRUE,
                                         ignore_pam=FALSE,
                                         cut_offset=NULL,
                                         standard_chr_only=TRUE,
                                         both_strands=TRUE,
                                         tss_window=NULL,
                                         n0_max=5,
                                         n1_max=100,
                                         n2_max=100
){
    tssObject   <- .validateTssObject(tssObject)
    guideSet    <- .validateGuideSet(guideSet)
    tss_window  <- .validateTssWindow(tss_window)
    seqlevelsStyle <- match.arg(seqlevelsStyle)
    guideSetTemp <- addSpacerAlignments(guideSet, 
                                        aligner=aligner,
                                        columnName=columnName,
                                        addSummary=TRUE,
                                        custom_seq=custom_seq,
                                        bowtie_index=bowtie_index,
                                        bwa_index=bwa_index,
                                        seqlevelsStyle=seqlevelsStyle,
                                        bsgenome=bsgenome,
                                        n_mismatches=0,
                                        all_alignments=FALSE, #To reduce burden
                                        canonical=canonical,
                                        ignore_pam=ignore_pam,
                                        cut_offset=cut_offset,
                                        standard_chr_only=standard_chr_only,
                                        both_strands=both_strands,
                                        tss_window=tss_window,
                                        txObject=txObject,
                                        tssObject=tssObject)

    cols <- .aln_cols(n_mismatches=n_mismatches)
            
    for (col in cols){
        mcols(guideSetTemp)[[col]] <- as.numeric(mcols(guideSetTemp)[[col]])
    }
    #Not calculating off-targets for guides with more than 5 on-targets:
    good <- which(guideSetTemp$n0<=n0_max)
    if (length(good)>0 & n_mismatches>=1){
        mcols(guideSetTemp)[good,] <- mcols(addSpacerAlignments(guideSetTemp[good],
                                                                aligner=aligner,
                                                                columnName=columnName,
                                                                addSummary=TRUE,
                                                                txObject=txObject,
                                                                tssObject=tssObject,
                                                                custom_seq=custom_seq,
                                                                bowtie_index=bowtie_index,
                                                                bwa_index=bwa_index,
                                                                seqlevelsStyle=seqlevelsStyle,
                                                                bsgenome=bsgenome, 
                                                                n_mismatches=1,
                                                                all_alignments=all_alignments,
                                                                canonical=canonical,
                                                                ignore_pam=ignore_pam,
                                                                cut_offset=cut_offset,
                                                                standard_chr_only=standard_chr_only,
                                                                both_strands=both_strands,
                                                                tss_window=tss_window))
    }
    good <- which(guideSetTemp$n0<=n0_max & guideSetTemp$n1<=n1_max)
    if (length(good)>0 & n_mismatches>=2){
        mcols(guideSetTemp)[good,] <- mcols(addSpacerAlignments(guideSetTemp[good],
                                                                aligner=aligner,
                                                                columnName=columnName,
                                                                addSummary=TRUE,
                                                                txObject=txObject,
                                                                tssObject=tssObject,
                                                                custom_seq=custom_seq,
                                                                bowtie_index=bowtie_index,
                                                                bwa_index=bwa_index,
                                                                seqlevelsStyle=seqlevelsStyle,
                                                                bsgenome=bsgenome, 
                                                                n_mismatches=2,
                                                                all_alignments=all_alignments,
                                                                canonical=canonical,
                                                                ignore_pam=ignore_pam,
                                                                cut_offset=cut_offset,
                                                                standard_chr_only=standard_chr_only,
                                                                both_strands=both_strands,
                                                                tss_window=tss_window))
    }
    good <- which(guideSetTemp$n0<=n0_max & guideSetTemp$n2<=n2_max)
    if (length(good)>0 & n_mismatches>=3){
        mcols(guideSetTemp)[good,] <- mcols(addSpacerAlignments(guideSetTemp[good],
                                                                aligner=aligner,
                                                                columnName=columnName,
                                                                addSummary=TRUE,
                                                                txObject=txObject,
                                                                tssObject=tssObject,
                                                                custom_seq=custom_seq,
                                                                bowtie_index=bowtie_index,
                                                                bwa_index=bwa_index,
                                                                seqlevelsStyle=seqlevelsStyle,
                                                                bsgenome=bsgenome, 
                                                                n_mismatches=3,
                                                                all_alignments=all_alignments,
                                                                canonical=canonical,
                                                                ignore_pam=ignore_pam,
                                                                cut_offset=cut_offset,
                                                                standard_chr_only=standard_chr_only,
                                                                both_strands=both_strands,
                                                                tss_window=tss_window))
    }
    mcols(guideSet)[[columnName]] <- mcols(guideSetTemp)[[columnName]]
    if (addSummary){
        mcols(guideSet)[cols] <- mcols(guideSetTemp)[cols]
    }
    return(guideSet)
}



.aln_cols <- function(n_mismatches,
                      excludePromoter=FALSE
){
    if (n_mismatches<=3){
        n_mismatches=3
    }
    cols <- c()
    for (k in c(0,seq_len(n_mismatches))){
        temp <- paste0("n", k, c("", "_c", "_p"))
        cols <- c(cols, temp)
    }
    
    if (excludePromoter){
        cols <- cols[!grepl("_p", cols)]
    }
    return(cols)
}





#' @rdname addSpacerAlignments
#' @export
#' @importFrom S4Vectors split
addSpacerAlignments <- function(guideSet,
                                aligner=c("bowtie", "bwa", "biostrings"),
                                columnName="alignments",
                                addSummary=TRUE,
                                txObject=NULL,
                                tssObject=NULL,
                                custom_seq=NULL,
                                bowtie_index=NULL,
                                bwa_index=NULL,
                                seqlevelsStyle=c("UCSC", "NCBI"),
                                bsgenome=NULL,
                                n_mismatches=0,
                                n_max_alignments=1000,
                                all_alignments=TRUE,
                                canonical=TRUE,
                                ignore_pam=FALSE,
                                cut_offset=NULL,
                                standard_chr_only=TRUE,
                                both_strands=TRUE,
                                tss_window=NULL
){
    aligner <- match.arg(aligner)
    tssObject <- .validateTssObject(tssObject)
    guideSet  <- .validateGuideSet(guideSet)
    crisprNuclease <- crisprNuclease(guideSet)
    seqlevelsStyle <- match.arg(seqlevelsStyle)
    tss_window   <- .validateTssWindow(tss_window)
    n_mismatches <- .validateNumberOfMismatches(n_mismatches,
                                                isBowtie=aligner=="bowtie")
    spacers <- spacers(guideSet,
                       as.character=TRUE)
    uniqueSpacers <- unique(spacers)
    aln <- getSpacerAlignments(uniqueSpacers,
                               aligner=aligner,
                               n_mismatches=n_mismatches,
                               custom_seq=custom_seq,
                               bowtie_index=bowtie_index,
                               bwa_index=bwa_index,
                               seqlevelsStyle=seqlevelsStyle,
                               bsgenome=bsgenome,
                               n_max_alignments=n_max_alignments,
                               all_alignments=all_alignments,
                               crisprNuclease=crisprNuclease,
                               canonical=canonical,
                               ignore_pam=ignore_pam,
                               cut_offset=cut_offset,
                               standard_chr_only=standard_chr_only,
                               both_strands=both_strands)
    metadata(aln)[["genome"]] <- metadata(guideSet)[["genome"]]


    if (!is.null(txObject)){
        aln <- .addGeneAnnotationColumns(aln,
                                         txObject=txObject)
    }
    if (!is.null(tssObject)){
        aln <- .addPromoterAnnotationColumns(aln,
                                             tssObject=tssObject,
                                             tss_window=tss_window)
    }
    if (addSummary){
        summary <- .getAlignmentsSummary(aln=aln,
                                         possibleGuides=spacers)           
        summary <- summary[, !grepl("_nc$", colnames(summary))]
        summary <- summary[match(spacers, summary$spacer),, drop=FALSE]
        cols <- .aln_cols(n_mismatches)
        if (sum(cols %in% colnames(mcols(guideSet)))>0){
            warning("Overwriting existing alignments summary. To ",
                    "avoid overwriting, set addSummary=FALSE.")
        }
        for (col in cols){
            mcols(guideSet)[[col]] <- summary[[col]]
        }

    }

    dfs <- S4Vectors::split(aln,
                            f=factor(aln$query,
                                     levels=uniqueSpacers))
    dfs <- dfs[spacers]
    names(dfs) <- names(guideSet)
    mcols(guideSet)[[columnName]] <- dfs
    return(guideSet)
}








#' @rdname addSpacerAlignments
#' @export
#' @importFrom GenomeInfoDb seqlevelsStyle<-
getSpacerAlignments <- function(spacers,
                                aligner=c("bowtie", "bwa", "biostrings"),
                                custom_seq=NULL,
                                bowtie_index=NULL,
                                bwa_index=NULL,
                                seqlevelsStyle=c("UCSC", "NCBI"),
                                bsgenome=NULL,
                                n_mismatches=0,
                                n_max_alignments=1000,
                                all_alignments=TRUE,
                                crisprNuclease=NULL,
                                canonical=TRUE,
                                ignore_pam=FALSE,
                                cut_offset=NULL,
                                standard_chr_only=TRUE,
                                both_strands=TRUE
){
    if (is(spacers, "GuideSet")){
        spacers <- spacers(spacers)
    }
    aligner <- match.arg(aligner)
    if (!is.null(custom_seq) & aligner != "biostrings"){
        aligner <- "biostrings"
        message("Setting aligner to biostrings since a custom_seq is provided.")

    }
    spacers <- as.character(spacers)
    seqlevelsStyle <- match.arg(seqlevelsStyle)
    if (aligner=="bowtie"){
        results <- .getSpacerAlignments_bowtie(spacers,
                                               bowtie_index=bowtie_index,
                                               bsgenome=bsgenome,
                                               n_mismatches=n_mismatches,
                                               n_max_alignments=n_max_alignments,
                                               all_alignments=all_alignments,
                                               crisprNuclease=crisprNuclease,
                                               canonical=canonical,
                                               ignore_pam=ignore_pam,
                                               cut_offset=cut_offset,
                                               standard_chr_only=standard_chr_only)
    } else if (aligner=="bwa"){
        results <- .getSpacerAlignments_bwa(spacers,
                                            bwa_index=bwa_index,
                                            bsgenome=bsgenome,
                                            n_mismatches=n_mismatches,
                                            crisprNuclease=crisprNuclease,
                                            canonical=canonical,
                                            ignore_pam=ignore_pam,
                                            cut_offset=cut_offset,
                                            standard_chr_only=standard_chr_only)
    } else if (aligner=="biostrings"){
        results <- .getSpacerAlignments_biostrings(spacers,
                                                   custom_seq=custom_seq,
                                                   n_mismatches=n_mismatches, 
                                                   crisprNuclease=crisprNuclease,
                                                   canonical=canonical,
                                                   both_strands=both_strands,
                                                   cut_offset=cut_offset)
    }
    if (aligner %in% c("bowtie", "bwa")){
        seqlevelsStyle(results) <- seqlevelsStyle
    }
    return(results)
}






#' @importFrom crisprBowtie runCrisprBowtie
#' @importFrom crisprBase pamLength pamSide
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb isCircular<- isCircular
#' @importFrom GenomeInfoDb seqinfo<- seqinfo Seqinfo
#' @importFrom GenomeInfoDb genome<- genome
.getSpacerAlignments_bowtie <- function(spacers,
                                        bowtie_index=NULL,
                                        bsgenome=NULL,
                                        n_mismatches=0,
                                        n_max_alignments=1000,
                                        all_alignments=TRUE,
                                        crisprNuclease=NULL,
                                        canonical=TRUE,
                                        ignore_pam=FALSE,
                                        cut_offset=NULL,
                                        standard_chr_only=TRUE
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    spacers    <- as.character(spacers)
    spacer_len <- unique(nchar(spacers))
    pamside <- pamSide(crisprNuclease)
    pamlen  <- pamLength(crisprNuclease)
    if (length(spacer_len) > 1){
        stop("All spacer sequences must have the same length.")
    }

    #Performing alignment:
    results <- runCrisprBowtie(spacers=spacers,
                               bowtie_index=bowtie_index,
                               n_mismatches=n_mismatches,
                               n_max_alignments=n_max_alignments,
                               crisprNuclease=crisprNuclease,
                               canonical=canonical,
                               ignore_pam=ignore_pam,
                               all_alignments=all_alignments,
                               force_spacer_length=TRUE)
    results$mmnuc1 <- NULL
    results$mmnuc2 <- NULL
    results$mmnuc3 <- NULL

    #If no alignments, return empty data frame:
    if (nrow(results)==0){
        emptyDF <- .returnEmptyAlignments(spacers,
                                          n_mismatches=n_mismatches)
        return(emptyDF)
    }

    #Filter and rename chromosomes:
    results <- results[results$chr!="chrEBV",,drop=FALSE]
    results$chr[results$chr=="chrMT"] <- "chrM"
    results$cut_site <- getCutSiteFromPamSite(pam_site=results$pam_site,
                                              strand=results$strand,
                                              crisprNuclease=crisprNuclease,
                                              cut_offset=cut_offset)
    #Setting mismatch locations relative to PAM site:
    wh <- which(results$n_mismatches!=0)
    if (pamside=="3prime"){
        rel.pam <- -(spacer_len + 1)
    } else {
        rel.pam <- pamlen-1
    }
    results$mm1[wh] <- results$mm1[wh] + rel.pam
    results$mm2[wh] <- results$mm2[wh] + rel.pam
    results$mm3[wh] <- results$mm3[wh] + rel.pam


    #Transforming to a GRanges:
    results <- .as_gr(results)
    results$query  <- DNAStringSet(results$spacer)
    results$spacer <- DNAStringSet(results$protospacer)
    results$protospacer <- NULL
    results$pam         <- DNAStringSet(results$pam)


    # Adding metadata:
    metadata(results)$crisprNuclease <- crisprNuclease
    metadata(results)$aligner  <- "bowtie"
    if (exists("n_mismatches")){
        metadata(results)$n_mismatches <- n_mismatches
    }
    if (exists("n_max_alignments")){
        metadata(results)$n_max_alignments <- n_max_alignments
    }
    if (exists("all_alignments")){
        metadata(results)$all_alignments <- all_alignments
    } 
    metadata(results)$spacer_len <- unique(nchar(as.character(results$spacer)))
    names(results) <- paste0("aln_", seq_along(results))



    if (!is.null(bsgenome)){
        info <- seqinfo(bsgenome)
    } else {
        levs <- unique(as.character(seqnames(results)))
        info <- GenomeInfoDb::Seqinfo(levs)
    }


    if (standard_chr_only){
        info <- keepStandardChromosomes(info)
        valid <- as.character(seqnames(results)) %in% seqnames(info)
        results <- results[valid]
        results <- keepStandardChromosomes(results)
    }
    
    if (!is.null(bsgenome)){
        genome(results)  <- genome(info)[1]
        lens  <- seqlengths(info)[as.character(seqnames(seqinfo(results)))]
        circs <- isCircular(info)[as.character(seqnames(seqinfo(results)))]
        seqlengths(seqinfo(results)) <- lens
        isCircular(seqinfo(results)) <- circs
    }
     
    return(results)
}




#' @importFrom crisprBwa runCrisprBwa
#' @importFrom crisprBase pamLength pamSide
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomeInfoDb isCircular<- isCircular
#' @importFrom GenomeInfoDb seqinfo<- seqinfo Seqinfo
#' @importFrom GenomeInfoDb genome<- genome
.getSpacerAlignments_bwa <- function(spacers,
                                     bwa_index=NULL,
                                     bsgenome=NULL,
                                     n_mismatches=0,
                                     crisprNuclease=NULL,
                                     canonical=TRUE,
                                     ignore_pam=FALSE,
                                     cut_offset=NULL,
                                     standard_chr_only=TRUE
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    spacers    <- as.character(spacers)
    spacer_len <- unique(nchar(spacers))
    pamside <- pamSide(crisprNuclease)
    pamlen  <- pamLength(crisprNuclease)
    if (length(spacer_len) > 1){
        stop("All spacer sequences must have the same length.")
    }

    #Performing alignment:
    results <- runCrisprBwa(spacers=spacers,
                            bwa_index=bwa_index,
                            n_mismatches=n_mismatches,
                            crisprNuclease=crisprNuclease,
                            canonical=canonical,
                            ignore_pam=ignore_pam,
                            force_spacer_length=TRUE)

    #If no alignments, return empty data frame:
    if (nrow(results)==0){
        emptyDF <- .returnEmptyAlignments(spacers,
                                          n_mismatches=n_mismatches)
        return(emptyDF)
    }

    #Filter and rename chromosomes:
    results <- results[results$chr!="chrEBV",,drop=FALSE]
    results$chr[results$chr=="chrMT"] <- "chrM"
    results$cut_site <- getCutSiteFromPamSite(pam_site=results$pam_site,
                                              strand=results$strand,
                                              crisprNuclease=crisprNuclease,
                                              cut_offset=cut_offset)
    
    #Setting mismatch locations relative to PAM site:
    wh <- which(results$n_mismatches!=0)
    if (pamside=="3prime"){
        rel.pam <- -(spacer_len + 1)
    } else {
        rel.pam <- pamlen-1
    }
    if (n_mismatches>0){
        for (k in seq_len(n_mismatches)){
            col <- paste0("mm", k)
            results[[col]][wh] <- results[[col]][wh] + rel.pam
        }
    }
  
    #Transforming to a GRanges:
    results <- .as_gr(results)
    results$query  <- DNAStringSet(results$spacer)
    results$spacer <- DNAStringSet(results$protospacer)
    results$protospacer <- NULL
    results$pam         <- DNAStringSet(results$pam)

    # Adding metadata:
    metadata(results)$crisprNuclease <- crisprNuclease
    metadata(results)$aligner  <- "bwa"
    if (exists("n_mismatches")){
        metadata(results)$n_mismatches <- n_mismatches
    } 
    metadata(results)$spacer_len <- unique(nchar(as.character(results$spacer)))
    names(results) <- paste0("aln_", seq_along(results))



    if (!is.null(bsgenome)){
        info <- seqinfo(bsgenome)
    } else {
        levs <- unique(as.character(seqnames(results)))
        info <- GenomeInfoDb::Seqinfo(levs)
    }


    if (standard_chr_only){
        info <- keepStandardChromosomes(info)
        valid <- as.character(seqnames(results)) %in% seqnames(info)
        results <- results[valid]
        results <- keepStandardChromosomes(results)
    }
    
    if (!is.null(bsgenome)){
        genome(results)  <- genome(info)[1]
        lens  <- seqlengths(info)[as.character(seqnames(seqinfo(results)))]
        circs <- isCircular(info)[as.character(seqnames(seqinfo(results)))]
        seqlengths(seqinfo(results)) <- lens
        isCircular(seqinfo(results)) <- circs
    }
     
    return(results)
}





#spacers = c("ACAAAACTGTGCTAGACATG","ACTAAACTGTGCTAGACAAC")
#custom_seq = "ACGAAACTCTGCTAGACATGTGGCGGTTTAGCCAGCTCCCCACATGTCTAGCACAGTTTTGTATGTAT"
#getSpacerAlignments(spacers, n_mismatches=5, nuclease="g", custom_seq=custom_seq, canonical=FALSE)
#getSpacerAlignments(spacers, n_mismatches=5, nuclease="Cas9", custom_seq=custom_seq, canonical=FALSE, both_strands=TRUE)
#spacers = c("ACAAAACTGTGCTAGACATGAGG","ACAAAACTATGCTAGACATGAGG", "ACAAAACTATAATAGACATGAGG")
#custom_seq = "ACGAATTTCACAAAACTGTGCTAGACATGAGGCTCCCCACATGTCCTCATGTCTAGCACAGTTTTGTCAAA"
#getSpacerAlignments(spacers, n_mismatches=5, nuclease="Cas12a", custom_seq=custom_seq, canonical=FALSE)
#getSpacerAlignments(spacers, n_mismatches=5, nuclease="Cas12a", custom_seq=custom_seq, canonical=FALSE, both_strands=TRUE)
#' @importFrom Biostrings matchPattern DNAStringSet
#' @importFrom IRanges ranges
#' @importFrom crisprBase pams pamLength pamSide
#' @importFrom utils adist
#' @importFrom GenomeInfoDb seqlengths<- isCircular<- genome<-
.getSpacerAlignments_biostrings <- function(spacers,
                                            custom_seq,
                                            n_mismatches=0, 
                                            crisprNuclease=NULL,
                                            canonical=TRUE,
                                            both_strands=TRUE,
                                            cut_offset=NULL
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    custom_seq <- .validateDNACharacterVariable(seq=custom_seq,
                                                argument="custom_seq",
                                                len=1,
                                                nullOk=FALSE,
                                                exactBases=FALSE)

    # Iterating over spacers:
    results <- lapply(spacers, function(spacer){
        locs_fwd <-  Biostrings::matchPattern(spacer,
                                              custom_seq,
                                              max.mismatch=n_mismatches)
        locs_fwd <- locs_fwd[BiocGenerics::start(locs_fwd)>0]
        locs_fwd <- locs_fwd[BiocGenerics::end(locs_fwd)<=nchar(custom_seq)]
        locs_fwd <- .view2df(locs_fwd,
                             crisprNuclease=crisprNuclease,
                             spacer=spacer,
                             strand="+",
                             custom_seq=custom_seq)

        locs <- locs_fwd
        if (both_strands){
            locs_rev <-  Biostrings::matchPattern(spacer,
                                                  .revComp(custom_seq),
                                                  max.mismatch=n_mismatches)
            locs_rev <- locs_rev[BiocGenerics::start(locs_rev)>0]
            locs_rev <- locs_rev[BiocGenerics::end(locs_rev)<=nchar(custom_seq)]
            locs_rev <- .view2df(locs_rev,
                                 crisprNuclease=crisprNuclease,
                                 spacer=spacer,
                                 strand="-",
                                 custom_seq=custom_seq)
            locs <- rbind(locs, locs_rev)
        }
        locs <- .removeInvalidProtospacers(df=locs,
                                           canonical=canonical,
                                           crisprNuclease=crisprNuclease,
                                           custom_seq=custom_seq,
                                           n_mismatches=n_mismatches)
        if (is.null(locs)){
            locs <- .returnEmptyAlignments(spacers,
                                           n_mismatches=n_mismatches)
        }
        return(locs)
    })

    # Combining
    ns <- vapply(results, nrow, FUN.VALUE=1)
    if (sum(ns==0)==length(ns)){
        locs <- .returnEmptyAlignments(spacers,
                                       n_mismatches=n_mismatches)
        return(locs)
    } else {
        locs <- results[ns>0]
        locs <- do.call(rbind, locs)
    }
    locs <- .annotateMismatches(locs,
                                n_mismatches=n_mismatches)

    #Transforming to a GRanges:
    locs <- .as_gr(locs)
    locs$query  <- DNAStringSet(locs$spacer)
    locs$spacer <- DNAStringSet(locs$protospacer)
    locs$protospacer <- NULL
    locs$pam         <- DNAStringSet(locs$pam)
    metadata(locs)$crisprNuclease   <- crisprNuclease
    metadata(locs)$genome     <- "custom"
    metadata(locs)$aligner    <- "biostrings"
    metadata(locs)$custom_seq <- DNAStringSet(custom_seq)
    metadata(locs)$spacer_len <- unique(nchar(as.character(locs$spacer)))
    if (exists("n_mismatches")) {
        metadata(locs)$n_mismatches <- n_mismatches
    }
    if (exists("both_strands")){
        metadata(locs)$both_strands <- both_strands  
    }
    names(locs) <- paste0("aln_", seq_along(locs))

    locs$cut_site <- getCutSiteFromPamSite(pam_site=locs$pam_site,
                                           strand=as.character(strand(locs)),
                                           crisprNuclease=crisprNuclease,
                                           cut_offset=cut_offset)
    #Taking care of seqinfo:
    lens <- width(metadata(locs)$custom_seq)
    names(lens) <- names(metadata(locs)$custom_seq)
    seqlengths(locs) <- lens
    isCircular(locs) <- FALSE
    genome(locs) <- "custom"
    return(locs)
}


  

.annotateMismatches <- function(df,
                                n_mismatches
){
    cols <- c("n_mismatches",
              paste0("mm", seq_len(n_mismatches)))
    df[cols] <- NA
    spacer_col <- "spacer"
    protospacer_col <- "protospacer"
    words1 <- df[[spacer_col]]
    words2 <- df[[protospacer_col]]
    df[["n_mismatches"]] <- vapply(seq_len(nrow(df)), function(i){
        utils::adist(words1[i], words2[i])[[1]]
    }, FUN.VALUE=1)
    words1 <- Biostrings::DNAStringSet(words1)
    words2 <- Biostrings::DNAStringSet(words2)
    x1 <- as.matrix(words1)
    x2 <- as.matrix(words2)
    whs <- apply(x1!=x2,1,which, simplify=FALSE)
    mm <- lapply(whs, function(wh){
        wh <- c(wh, rep(NA, n_mismatches-length(wh)))
        return(wh)
    }) 
    mm <- do.call(rbind, mm)
    df[, paste0("mm", seq_len(n_mismatches))] <- mm
    return(df)
}


# Remove targets with invalid PAM sequences,
# or outside of custom sequence. 
.removeInvalidProtospacers <- function(df,
                                       canonical,
                                       custom_seq,
                                       crisprNuclease,
                                       n_mismatches
){
    if (is.null(df)) {
        return(df)
    }
    pamlen <- pamLength(crisprNuclease)
    seqlen <- nchar(custom_seq)
    pamside <- pamSide(crisprNuclease)
    if (pamside=="3prime"){
        good_fwd <- df$strand=="+" & (df$pam_site + pamlen-1 <=seqlen)
        good_rev <- df$strand=="-" & (df$pam_site - pamlen+1 >=1)
    } else {
        good_fwd <- df$strand=="+" & (df$pam_site >= 1)
        good_rev <- df$strand=="-" & (df$pam_site <= seqlen)
    }
    good <- good_fwd | good_rev
    if (sum(good)==0){
        df <- .returnEmptyAlignments(spacers,
                                     n_mismatches=n_mismatches)
        return(df)
    } else {
        df <- df[good,,drop=FALSE]
    }
    df$pam <- getPAMSequence_customSeq(custom_seq=custom_seq,
                                       pam_site=df$pam_site,
                                       strand=df$strand,
                                       crisprNuclease=crisprNuclease)
    df$canonical <- df$pam %in% pams(crisprNuclease,
                                     primary=TRUE,
                                     as.character=TRUE)
    if (canonical){
        df <- df[df$canonical,,drop=FALSE]
        if (nrow(df)<1){
            df <- NULL
        }
    }
    return(df)
}



# Convert view object to a data.frame
.view2df <- function(view,
                     strand=c("+", "-"),
                     crisprNuclease,
                     spacer,
                     custom_seq
){
    seqlen <- nchar(custom_seq)
    strand <- match.arg(strand)
    if (length(view)==0){
        df <- NULL
        return(df)
    } 

    df <- as.data.frame(IRanges::ranges(view))
    df$strand <- strand
    seqName <- names(custom_seq)
    if (is.null(seqName)){
        seqName <- "customSequence"
    }
    df$chr <- seqName
    spacer_len <- unique(df$width)
    if (strand=="-"){
        start <- seqlen-df$end+1
        end   <- seqlen-df$start+1
        df$start <- start
        df$end   <- end
    }
    df$pam_site <- getPAMSiteFromStartAndEnd(start=df$start,
                                             end=df$end,
                                             crisprNuclease=crisprNuclease,
                                             strand=df$strand)

    df$start <- df$end <- df$width <- NULL
    df[["protospacer"]] <- as.character(view)
    df[["spacer"]] <- spacer
    cols <- c("chr", "pam_site", "strand")
    cols <- c(cols, setdiff(colnames(df), cols))
    df <- df[,cols,drop=FALSE]
    return(df)
}








.returnEmptyAlignments <- function(spacers,
                                   n_mismatches
){
    spacer_len <- unique(nchar(spacers))
    spacer_col <- paste0('spacer_', spacer_len, 'mer')
    protospacer_col <- paste0('protospacer_', spacer_len, 'mer')
    emptyDF <- data.frame(chr=character(0),
                          strand=character(0),
                          pam_site=numeric(0),
                          pam=character(0),
                          cut_site=numeric(0),
                          canonical=logical(0),
                          n_mismatches=numeric(0))
    if (n_mismatches>0){
        for (k in seq_len(n_mismatches)){
            col <- paste0("mm", k)
            emptyDF[[col]] <- numeric(0)
        }
    }

    emptyDF[[spacer_col]] <- character(0)
    emptyDF[[protospacer_col]] <- character(0)
    cols <- c((ncol(emptyDF)-1), ncol(emptyDF), seq_len(ncol(emptyDF)-2))
    emptyDF <- emptyDF[,cols]
    return(emptyDF)
}




.getAlignmentsSummary <- function(aln,
                                  possibleGuides=NULL
){

    max_mm=metadata(aln)$n_mismatches
    if (.isGRanges(aln)){
        aln <- .as_df(aln)
    }
    stopifnot(is.data.frame(aln))

    if (nrow(aln)==0){
        if (is.null(possibleGuides)){
            stop('No alignments to summarize. Pass input guides",
                 " to "possibleGuides" to avoid this error.')
        }
        warning('No alignments found, still adding summary data.')
        # return "empty" alignment summary (all values set to 0)
        df <- data.frame(spacer=possibleGuides)
        cols <- c('nx', 'nx_c', 'nx_nc', 'nx_p')
        for (i in cols){
            for (j in 0:max_mm){
                df[[gsub('x', j, i)]] <- 0
            }
        }
        return(df)
    }

    df <- data.frame(spacer=unique(aln$query))

    # Create summary columns:
    cols <- c('', '_c', '_nc', '_p')
    nn <- ifelse(max_mm<=3,3, max_mm)
    cols <- lapply(cols, function(x){
        paste0('n', c(0,seq_len(nn)), x)
    })
    cols <- unlist(cols)
    df[cols] <- NA
    default.p <- ifelse('promoters' %in% colnames(aln), 0, NA)
    df[grep(paste0('^n[0-', max_mm, ']$'), cols, value=TRUE)] <- 0
    df[grep(paste0('^n[0-', max_mm, ']\\_n?c$'), cols, value=TRUE)] <- 0
    df[grep(paste0('^n[0-', max_mm, ']\\_p$'), cols, value=TRUE)]  <- default.p

    # return default values if there are no alignments
    if (sum(aln$pam_site==0)==nrow(aln)){
        return(df)
    }

    # count number of alignments per column
    .countAlignments <- function(col, cond){
        mm <- as.integer(substr(col, start=2, stop=2))
        tab <- table(aln$query[aln$n_mismatches==mm & cond])
        if (length(tab) > 0){
            wh <- match(names(tab), df$spacer)
            temp <- df[,col]
            temp[wh] <- tab
            return(temp)
        } else {
            return(df[,col])
        }
    }

    # loop over columns
    for (col in colnames(df)[-1]){
        if (sum(is.na(df[[col]]))==0){
            aln.type <- gsub('[0-9]', '', col)
            cond <- switch(EXPR=aln.type,
                           'n'=TRUE,
                           'n_c'=!is.na(aln$cds),
                           'n_nc'=is.na(aln$cds),
                           'n_p'=!is.na(aln$promoters))
            df[[col]] <- .countAlignments(col=col,
                                          cond=cond)
        }
    }

    # add non-aligning guides, if any
    if (!is.null(possibleGuides)){
        missing <- setdiff(possibleGuides, df$spacer)
        if (length(missing)>0){
            temp <- data.frame(spacer=missing)
            cols <- setdiff(colnames(df), "spacer")
            mat <- matrix(0, length(missing), length(cols))
            temp <- cbind(temp, mat)
            colnames(temp)[-1] <- cols
            df <- rbind(df, temp)
        }
    }
    return(df)
}





.validateNumberOfMismatches <- function(n_mismatches,
                                        isBowtie=FALSE
){
    stopifnot(is.numeric(n_mismatches))
    if (isBowtie){
        if (!n_mismatches %in% c(0, 1, 2, 3)){
            stop('For bowtie, "n_mismatches" must be an integer ',
                 'between 0 and 3, inclusive.')
        }
    }
    return(n_mismatches)
}











## This function adds gene annotation columns to the alignments 
## obtained from getSpacerAlignments.
##
#' @importFrom S4Vectors metadata mcols mcols<- subjectHits queryHits
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom IRanges IRanges ranges ranges<- findOverlaps
#' @importFrom IRanges precede follow 
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
.addGeneAnnotationColumns <- function(aln,
                                      anchor=c("cut_site", "pam_site"),
                                      txObject=NULL,
                                      ignore.strand=TRUE
){
    genome <- metadata(aln)$genome
    anchor <- match.arg(anchor)
    anchor <- .validateAnchor(anchor, aln)
    
    #To change...
    if (length(aln) == 0) {
        aln$cds         <- character(0)
        aln$fiveUTRs    <- character(0)
        aln$threeUTRs   <- character(0)
        aln$exons       <- character(0)
        aln$introns     <- character(0)
        aln$intergenics <- character(0)
        return(aln)
    }
    
    aln_gr <- aln
    ranges(aln_gr) <- IRanges(start=mcols(aln_gr)[[anchor]], width=1)
    ref <- .validateGRangesList(txObject)
    GenomeInfoDb::seqlevelsStyle(ref) <- GenomeInfoDb::seqlevelsStyle(aln)
    
    # Adding gene info columns:
    regions <- c("cds",
                 "fiveUTRs",
                 "threeUTRs",
                 "exons",
                 "introns")
    for (region in regions) {
        overlaps <- suppressWarnings(findOverlaps(aln_gr,
                                                  ref[[region]],
                                                  ignore.strand=ignore.strand))
        aln <- .addOverlapGenes(dat=aln,
                                ref=ref,
                                overlaps=overlaps,
                                region=region)
    }
    
    # Adding intergenic info when everything is empty:
    mcols(aln)$intergenic <- NA
    geneAnn <- mcols(aln[, regions, drop=FALSE])
    wh <- which(rowSums(is.na(geneAnn)) == ncol(geneAnn))
    if (length(wh) > 0) {
        upstream   <- precede(aln_gr[wh], ref$transcripts)
        downstream <- follow(aln_gr[wh], ref$transcripts)
        genes_up   <- ref$transcripts$gene_symbol[upstream]
        genes_up[is.na(genes_up)] <- ref$transcripts$tx_id[upstream][is.na(genes_up)]
        genes_down <- ref$transcripts$gene_symbol[downstream]
        genes_down[is.na(genes_down)] <- ref$transcripts$tx_id[downstream][is.na(genes_down)]
        flanking_genes <- paste(genes_up, "<->", genes_down)
        flanking_genes[flanking_genes == "NA <-> NA"] <- NA
        mcols(aln)$intergenic[wh] <- flanking_genes
    }
    return(aln)
}



.addOverlapGenes <- function(dat,
                             ref,
                             overlaps,
                             region
){
    mcols(dat)[[region]] <- NA
    if (length(overlaps) == 0) {
        return(dat)
    }
    if (.isGRangesList(ref)) {
        tx <- mcols(ref$transcripts)
        gr_gene <- mcols(ref[[region]])
    } else {
        tx <- NULL
        gr_gene <- mcols(ref)
    }
    if ("gene_symbol" %in% names(gr_gene)){
        gene_col <- "gene_symbol"
    } else if ("gene_id" %in% names(gr_gene)){
        gene_col <- "gene_id"
    } else {
        gene_col <- "tx_id"
    }
    mcols(overlaps)$gene <- gr_gene[[gene_col]][subjectHits(overlaps)]
    if (!is.null(tx)) {
        mcols(overlaps)$gene <- tx$gene_symbol[match(mcols(overlaps)$gene, 
                                                     tx[[gene_col]])]
    }
    genes <- split(mcols(overlaps)$gene,
                   f=queryHits(overlaps))
    genes <- vapply(genes, function(x) {
        paste0(unique(x), collapse = ";")
    }, FUN.VALUE="a")
    mcols(dat)[[region]][match(names(genes), seq_along(dat))] <- genes
    return(dat)
}






.addPromoterAnnotationColumns <- function(aln,
                                          anchor=c("cut_site","pam_site"),
                                          tssObject=NULL,
                                          tss_window=NULL,
                                          ignore.strand=TRUE
){
    tssObject <- .validateTssObject(tssObject)
    # genome <- metadata(aln)$genome
    anchor <- match.arg(anchor)
    anchor <- .validateAnchor(anchor, aln)
    tss_window <- .validateTssWindow(tss_window)
    
    #To change...
    if (length(aln) == 0) {
        aln$promoters <- character(0)
        return(aln)
    }
    
    aln_gr <- aln
    ranges(aln_gr) <- IRanges(start=mcols(aln_gr)[[anchor]], width=1)
    if (!is.null(tssObject)){
        ref_tss <- promoters(tssObject,
                             upstream = (-1 * tss_window[1]), 
                             downstream = tss_window[2])
    } else {
        ref_tss <- NULL
    }
    
    # Adding TSS info columns:
    if (!is.null(ref_tss)) {
        overlaps <- suppressWarnings(findOverlaps(aln_gr,
                                                  ref_tss,
                                                  ignore.strand=ignore.strand))
        
        aln <- .addOverlapGenes(dat=aln,
                                ref=ref_tss,
                                overlaps=overlaps,
                                region="promoters")
    }
    return(aln)
}
