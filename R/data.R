#' Example of a \linkS4class{TxDb} object converted to a
#'     \linkS4class{GRangesList}
#'
#' Example of a \linkS4class{TxDb} object converted to a
#'     \linkS4class{GRangesList} object for human gene
#'     IQSEC3 (ENSG00000120645).
#' 
#' @format Named \linkS4class{GRangesList} with 7 elements:
#'     \code{transcripts}, \code{exons},
#'     \code{cds}, \code{fiveUTRs}, \code{threeUTRs}, \code{introns}
#'     and \code{tss}.
#' 
#' @details The full human transcriptome \linkS4class{TxDb} object was
#'     obtained from the Ensembl 104 release using the \code{\link{getTxDb}}
#'     function and converted to a \linkS4class{GRangesList} object using the 
#'     \code{\link{TxDb2GRangesList}} function and subsequently subsetted
#'     to only contain the IQSEC3 gene (ENSG00000120645) located
#'     at the start of chr12 in the human genome (hg38 build).
#' @usage data(grListExample, package="crisprDesign")
"grListExample"





#' Example of a \linkS4class{GRanges} object containing TSS coordinates
#'
#' Example of a \linkS4class{GRanges} containing transcription starting
#'     site (TSS) coordinates for human gene IQSEC3 (ENSG00000120645).
#' 
#' @format \linkS4class{GRanges} object of length 2 corresponding to the 2 
#'    TSSs of gene IQSEC3.
#' 
#' @details The TSS coordinates were obtained from the two transcript
#'     stored in the \code{grListExample} object for gene IQSEC3.
#' @usage data(tssObjectExample, package="crisprDesign")
"tssObjectExample"



#' Example of a \linkS4class{GuideSet} object storing gRNA sequences targeting
#'     the CDS of IQSEC3
#'
#' Example of a \linkS4class{GuideSet} object storing gRNA sequences targeting
#'    the coding sequence of human gene IQSEC3 (ENSG00000120645) for SpCas9
#'    nuclease. 
#' 
#' @format A \linkS4class{GuideSet} object.
#' 
#' @details The object was obtained by calling \code{\link{findSpacers}} on the
#'     CDS region of human gene IQSEC3. See code in
#'     \code{inst/scripts/generateGuideSet.R}.
#' @usage data(guideSetExample, package="crisprDesign")
"guideSetExample"




#' Example of a fully-annotated \linkS4class{GuideSet} object storing gRNA
#'     sequences targeting the CDS of IQSEC3
#'
#' Example of a fully-annotated \linkS4class{GuideSet} object storing gRNA
#'    sequences targeting the coding sequence of human gene IQSEC3
#'    (ENSG00000120645) for SpCas9 nuclease. 
#' 
#' @format A \linkS4class{GuideSet} object.
#' 
#' @details The object was obtained by applying all available \code{add*}
#'     annotation functions in \code{crisprDesign} (e.g.
#'     \code{addSequenceFeatures}) to a randomly selected 20-guide subset of
#'      \code{guideSetExample}. See code in
#'     \code{inst/scripts/generateGuideSetFullAnnotation.R}.
#' @usage data(guideSetExampleFullAnnotation, package="crisprDesign")
"guideSetExampleFullAnnotation"




#' Example of a \linkS4class{GRanges} object containing repeat elements
#'
#' Example of a \linkS4class{GRanges} object containing genomic coordinates
#'     of repeat elements found in the neighborhood of human gene IQSEC3
#'     (ENSG00000120645).
#' 
#' @format A \linkS4class{GRanges} object. 
#' @usage data(grRepeatsExample, package="crisprDesign")
"grRepeatsExample"
