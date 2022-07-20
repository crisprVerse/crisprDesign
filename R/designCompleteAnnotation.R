# library(crisprDesign)
# library(crisprDesignGne)
# library(crisprBase)
# library(crisprDesignData)
# species="human"
# bowtie_index <- getBowtieIndex(species=species)
# bsgenome  <- getGenomePackage(species=species)
# snpFile <- getSNPFile() 
# data(SpCas9, package="crisprBase")
# data(CasRx, package="crisprBase")
# data(txdb_human, package="crisprDesignData")
# data(tss_human, package="crisprDesignData")
# data(mrnasHuman, package="crisprDesignData")
# mrnas <- mrnasHuman
# txObject = txdb_human
# queryValue <- c("ENST00000538872")
# bowtie_index <- getBowtieIndex(species=species, what="rna")


# gs <- precomputeGuides("ENSG00000133703",
#                        modality="crisprko",
#                        bsgenome=bsgenome,
#                        bowtie_index=bowtie_index,
#                        crisprNuclease=SpCas9,
#                        txObject=txdb_human,
#                        tssObject=tssdb_human,
#                        grRepeats=grRepeats)

# chromatinFiles <- getChromatinFiles()
# fastaFile <- getGenomeFasta()

# modality="CRISPRkd"
# crisprNuclease=CasRx
# canonical_ontarget=TRUE
# canonical_offtarget=FALSE
# all_alignments=TRUE 
# n_mismatches=1



#' @title One-step gRNA design and annotation function
#' 
#' @description One-step gRNA design and annotation function
#'    to faciliate the design and generation of genome-wide
#'    gRNA databases for a combination of parameters such 
#'    as nuclease, organism, and CRISPR modality.
#' 
#' @param queryValue Vector specifying the value(s) to search for in
#'     \code{txObject[[featureType]][[queryColumn]]}.
#' @param queryColumn Character string specifying the column in 
#'     \code{txObject[[featureType]]} to search for \code{queryValue}(s).
#' @param featureType For CRISPRko, string specifying the type of genomic
#'     feature to use to design gRNAs. Must be of the following:
#'     "transcripts", "exons", "cds", "fiveUTRs", "threeUTRs" or "introns".
#'     The default is "cds".
#' @param modality String specifying the CRISPR modality. Must be one of
#'     the following: "CRISPRko", "CRISPRa", "CRISPRi" or "CRISPRkd".
#'     CRISPRkd is reserved for DNA-targeting nucleases only such as CasRx.
#' @param bsgenome A \linkS4class{BSgenome} object from which to extract
#'     sequences if a \linkS4class{GRanges} object is provided as input. 
#' @param bowtie_index String specifying path to a bowtie index.
#' @param vcf Either a character string specfying a path to a VCF file
#'     or connection, or a \linkS4class{VCF} object.
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}} for annotating
#'     on-target and off-target alignments using gene annotation.
#' @param tssObject A \linkS4class{GRanges} object specifying TSS coordinates.
#' @param grRepeats A \linkS4class{GRanges} object containing repeat
#'     elements regions.
#' @param scoring_methods Character vector to specify which on-target scoring
#'     methods should be calculated. See crisprScore package to obtain
#'     available methods. 
#' @param tss_window Vector of length 2 specifying the start and coordinates
#'     of the CRISPRa/CRISPRi target region with respect to the TSS position.
#' @param n_mismatches Maximum number of mismatches permitted between guide RNA
#'     and genomic DNA.
#' @param max_mm The maximimum number of mismatches between a spacer and
#'     an off-target to be accepted when calculating aggregate off-target
#'     scores. 2 by default. 
#' @param canonical_ontarget Should only canonical PAM sequences be searched
#'     for designing gRNAs? TRUE by default.
#' @param canonical_offtarget Should only canonical PAM sequences by searched
#'     during the off-target search? TRUE by default.
#' @param all_alignments Should all all possible alignments be returned?
#'     TRUE by default.
#' @param fastaFile String specifying fasta file of the hg38 genome. Only 
#'     used for CRISPRa/i modality with hg38 genome and SpCas9 nuclease.
#'     This is needed to generate the CRISPRai scores. See the function
#'     \code{addCrispraiScores} for more details. 
#' @param chromatinFiles Named character vector of length 3 specifying
#'     BigWig files containing chromatin accessibility data. Only 
#'     used for CRISPRa/i modality with hg38 genome and SpCas9 nuclease.
#'     This is needed to generate the CRISPRai scores. See the function
#'     \code{addCrispraiScores} for more details. 
#' @param verbose Should messages be printed?
#' 
#' @return A \code{GuideSet} object.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom crisprBase spacerLength
#' @importFrom utils data
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @export
designCompleteAnnotation <- function(queryValue=NULL,
                                     queryColumn="gene_id",
                                     featureType="cds",
                                     modality=c("CRISPRko",
                                                "CRISPRa",
                                                "CRISPRi",
                                                "CRISPRkd"),
                                     bsgenome=NULL,
                                     bowtie_index=NULL,
                                     vcf=NULL,
                                     crisprNuclease=NULL,
                                     tssObject=NULL,
                                     txObject=NULL,
                                     grRepeats=NULL,
                                     scoring_methods=NULL,
                                     tss_window=NULL,
                                     n_mismatches=3,
                                     max_mm=2,
                                     canonical_ontarget=TRUE,
                                     canonical_offtarget=FALSE,
                                     all_alignments=TRUE,
                                     fastaFile=NULL,
                                     chromatinFiles=NULL,
                                     verbose=TRUE
){
    modality <- match.arg(modality)
    
    # Checking nucleases:
    data(SpCas9,     package="crisprBase", envir=environment())
    data(enAsCas12a, package="crisprBase", envir=environment())
    data(CasRx,      package="crisprBase", envir=environment())
    isCas9   <- .identicalNucleases(crisprNuclease, SpCas9)
    isCas12a <- .identicalNucleases(crisprNuclease, enAsCas12a)
    isCas13d <- .identicalNucleases(crisprNuclease, CasRx)
    isKD <- modality=="CRISPRkd"
    isKO <- modality=="CRISPRko"
    isA  <- modality=="CRISPRa"
    isI  <- modality=="CRISPRi"



    if (isA | isI){
        if (is.null(tss_window)){
            stop("tss_window must be specified for CRISPRa/i.")
        }
    }
    if (sum(c(isKD, isCas13d))==1){
        stop("Modality must be CRISPRkd when nuclease is set to CasRx, and vice versa.")        
    }

    if (isKO){
        gr <- queryTxObject(txObject,
                            featureType=featureType,
                            queryValue=queryValue,
                            queryColumn=queryColumn)
    } else if (isA | isI){
        gr <- queryTss(tssObject,
                       queryValue=queryValue,
                       queryColumn=queryColumn,
                       tss_window=tss_window)
    } else if (isKD){
        mrna <- getMrnaSequences(txids=queryValue,
                                 bsgenome=bsgenome,
                                 txObject=txObject)
        gr <- mrna
    }
    if (is.null(gr)){
        out <- NA
        return(out)
    }

    out <- findSpacers(x=gr,
                       crisprNuclease=crisprNuclease,
                       canonical=canonical_ontarget,
                       strict_overlap=FALSE,
                       bsgenome=bsgenome)
    out <- unique(out)
     # Renaming spacers:
    if (length(out)>0){
        names(out) <- paste0(queryValue, "_", seq_along(out))
    }
    if (is.null(out)){
        out <- NA
        return(out)
    } else if (!is.null(grRepeats) & !isKD){
        out <- removeRepeats(out,
                             gr.repeats=grRepeats)
        if (length(out)==0){
            out <- NULL
        }
    }
    # Cases where no guides are found:
    if (is.null(out)){
        out <- NA
        return(out)
    }

    # Cases where no guides are found in CDS:
    if (!isKD){
        chr <- as.character(seqnames(out))
        if (isKO){
            gr_grna <- GRanges(chr,IRanges(start=out$cut_site,
                                           end=out$cut_site))
        } else if (isA | isI){
            gr_grna  <- GRanges(chr, IRanges(start=out$pam_site,
                                             end=out$pam_site))
        } 
        overlaps <- findOverlaps(gr_grna, gr)
        if (length(overlaps)==0){
            out <- NA
            return(out)
        }
    }

    if (verbose){
        cat("[designCompleteAnnotation] Adding sequence statistics \n")
    }
    out <- addSequenceFeatures(out, addHairpin=FALSE)
    if (verbose){
        cat("[designCompleteAnnotation] Adding spacer alignments \n")
    }
    spacer_len <- spacerLength(crisprNuclease)
    good <- nchar(spacers(out, as.character=TRUE))==spacer_len
    out <- out[good]
    if (sum(good)==0){
        out <- NA
        return(out)
    }


 
    if (!isKD){
        out <- addSpacerAlignmentsIterative(out,
                                            aligner="bowtie",
                                            aligner_index=bowtie_index,
                                            n_mismatches=n_mismatches,
                                            canonical=canonical_offtarget,
                                            all_alignments=all_alignments,
                                            bsgenome=bsgenome,
                                            txObject=txObject,
                                            tssObject=tssObject)
    } else {
        out <- addSpacerAlignments(out,
                                   aligner="bowtie",
                                   aligner_index=bowtie_index,
                                   n_mismatches=n_mismatches,
                                   canonical=canonical_offtarget,
                                   all_alignments=all_alignments,
                                   bsgenome=NULL,
                                   txObject=txObject,
                                   tssObject=tssObject)
    }
   
    if (verbose){
        cat("[designCompleteAnnotation] Adding gene annotation \n")
    } 
    out <- addGeneAnnotation(out,
                             txObject=txObject)
    out <- addRestrictionEnzymes(out)

    if (isA | isI){
        if (verbose){
            cat("[designCompleteAnnotation] Adding TSS annotation \n")
        } 
        out <- addTssAnnotation(out,
                                tssObject=tssObject,
                                tss_window=tss_window)
    }

    if (verbose){
        cat("[designCompleteAnnotation] Adding on-target scores \n")
    } 
    out <- addOnTargetScores(out,
                             methods=scoring_methods)


    if (!is.null(fastaFile) & !is.null(chromatinFiles) & modality %in% c("CRISPRi", "CRISPRi") & isCas9){
        if (verbose){
            cat("[designCompleteAnnotation] Adding CRISPRai scores \n")
        }
        out <- addCrispraiScores(out,
                                 gr=gr,
                                 tssObject=tssObject,
                                 modality=modality,
                                 chromatinFiles=chromatinFiles,
                                 fastaFile=fastaFile)
    }

    if (isKO & isCas9){
        if (verbose){
            cat("[designCompleteAnnotation] Adding CFD scores annotation \n")
            out <- addOffTargetScores(out, max_mm=max_mm)
        }
    }

    if (!is.null(vcf) & !isKD){
        if (verbose){
            cat("[designCompleteAnnotation] Adding SNP annotation \n")
        }
        out <- addSNPAnnotation(out, vcf=vcf)
    }
    return(out)
}


