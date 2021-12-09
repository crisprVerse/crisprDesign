#' @param guideSet A \linkS4class{GuideSet} object.
#' @param aligner Which genomic alignment method should be used?
#'     Must be one of "bowtie", "blast", and "biostrings".
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
#' @param bowtie_index Path to the bowtie index to be used for alignment.
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



# addSpacerAlignmentsIterative <- function(guideSet,
#                                          aligner=c("bowtie", "blast", "biostrings"),
#                                          columnName="alignments",
#                                          addSummary=TRUE,
#                                          txObject=NULL,
#                                          tssObject=NULL,
#                                          custom_seq=NULL,
#                                          bowtie_index=NULL,
#                                          seqlevelsStyle=c("UCSC", "NCBI"),
#                                          bsgenome=NULL,
#                                          n_mismatches=0,
#                                          all_alignments=FALSE,
#                                          canonical=TRUE,
#                                          ignore_pam=FALSE,
#                                          cut_offset=NULL,
#                                          standard_chr_only=TRUE,
#                                          both_strands=TRUE,
#                                          tss_window=NULL,
#                                          n0_max=5,
#                                          n1_max=100,
#                                          n2_max=100


# many arguments to test













# addSpacerAlignments <- function(guideSet,
#                                 aligner=c("bowtie", "blast", "biostrings"),
#                                 columnName="alignments",
#                                 addSummary=TRUE,
#                                 txObject=NULL,
#                                 tssObject=NULL,
#                                 custom_seq=NULL,
#                                 bowtie_index=NULL,
#                                 seqlevelsStyle=c("UCSC", "NCBI"),
#                                 bsgenome=NULL,
#                                 n_mismatches=0,
#                                 n_max_alignments=1000,
#                                 all_alignments=TRUE,
#                                 canonical=TRUE,
#                                 ignore_pam=FALSE,
#                                 cut_offset=NULL,
#                                 standard_chr_only=TRUE,
#                                 both_strands=TRUE,
#                                 tss_window=NULL


# many args to test...














# getSpacerAlignments <- function(spacers,
#                                 aligner=c("bowtie", "blast", "biostrings"),
#                                 custom_seq=NULL,
#                                 bowtie_index=NULL,
#                                 seqlevelsStyle=c("UCSC", "NCBI"),
#                                 bsgenome=NULL,
#                                 n_mismatches=0,
#                                 n_max_alignments=1000,
#                                 all_alignments=TRUE,
#                                 crisprNuclease=NULL,
#                                 canonical=TRUE,
#                                 ignore_pam=FALSE,
#                                 cut_offset=NULL,
#                                 standard_chr_only=TRUE,
#                                 both_strands=TRUE



# many args to test...
