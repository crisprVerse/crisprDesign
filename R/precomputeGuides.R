

# bowtie_index <- getBowtieIndex(species=species)
# bsgenome  <- getGenomePackage(species=species)
# snpFile <- getSNPFile() 
# library(crisprBase)
# library(crisprDesignData)
# data(SpCas9, package="crisprBase")
# data(txdb_human, package="crisprDesignData")
# data(tss_human, package="crisprDesignData")
# repeats <- readRDS(file.path(geneDir, "repeats.rds"))
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




    




#' @importFrom crisprBase spacerLength
#' @importFrom utils data
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @export
precomputeGuides <- function(geneid,
                             modality=c("CRISPRko","CRISPRa", "CRISPRi"),
                             bsgenome=NULL,
                             bowtie_index=NULL,
                             vcf=NULL,
                             crisprNuclease=NULL,
                             tssObject=NULL,
                             txObject=NULL,
                             grRepeats=NULL,
                             scoring_methods=NULL,
                             verbose=TRUE,
                             tss_window=NULL,
                             n_mismatches=3,
                             max_mm=2,
                             canonical_ontarget=TRUE,
                             canonical_offtarget=FALSE,
                             all_alignments=TRUE,
                             remove_repeats=TRUE,
                             fastaFile=NULL,
                             chromatinFiles=NULL
){
    library(crisprScore)
    modality <- match.arg(modality)
    
    if (modality!="CRISPRko"){
        if (is.null(tss_window)){
            stop("tss_window must be specified for CRISPRa/i.")
        }
    }
    

    if (modality=='CRISPRko'){
        gr <- queryTxObject(txObject,
                            featureType="cds",
                            queryValue=geneid,
                            queryColumn="gene_id")
    } else if (modality=='CRISPRa' | modality=='CRISPRi'){
        gr <- queryTss(tssObject,
                       queryValue=geneid,
                       queryColumn="gene_id",
                       tss_window=tss_window)
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
        names(out) <- paste0(geneid, "_", seq_along(out))
    }
    if (is.null(out)){
        out <- NA
        return(out)
    } else if (remove_repeats){
        out <- removeRepeats(out, gr.repeats=grRepeats)
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
    chr <- as.character(seqnames(out))
    if (modality=="CRISPRko"){
        gr_grna <- GRanges(chr,IRanges(start=out$cut_site,
                                       end=out$cut_site))
    } else if (modality=="CRISPRa" | modality=="CRISPRi"){
        gr_grna  <- GRanges(chr, IRanges(start=out$pam_site,
                                         end=out$pam_site))
    }
    overlaps <- findOverlaps(gr_grna, gr)
    if (length(overlaps)==0){
        out <- NA
        return(out)
    }


    if (verbose){
        cat("[precomputeGuides] Adding sequence statistics \n")
    }
    out <- addSequenceFeatures(out, addHairpin=FALSE)
    if (verbose){
        cat("[precomputeGuides] Adding spacer alignments \n")
    }
    spacer_len <- spacerLength(crisprNuclease)
    good <- nchar(spacers(out, as.character=TRUE))==spacer_len
    out <- out[good]
    if (sum(good)==0){
        out <- NA
        return(out)
    }
    out <- addSpacerAlignmentsIterative(out,
                                        aligner="bowtie",
                                        aligner_index=bowtie_index,
                                        n_mismatches=n_mismatches,
                                        canonical=canonical_offtarget,
                                        all_alignments=all_alignments,
                                        bsgenome=bsgenome,
                                        txObject=txObject,
                                        tssObject=tssObject)
    if (verbose){
        cat("[precomputeGuides] Adding gene annotation \n")
    } 
    out <- addGeneAnnotation(out, 
                             txObject=txObject)
    out <- addRestrictionEnzymes(out)

    if (modality=='CRISPRa' | modality=='CRISPRi'){
        if (verbose){
            cat("[precomputeGuides] Adding TSS annotation \n")
        } 
        out <- addTssAnnotation(out,
                                tssObject=tssObject,
                                tss_window=tss_window)
    }

    if (verbose){
        cat("[precomputeGuides] Adding on-target scores \n")
    } 
    out <- addOnTargetScores(out,
                             methods=scoring_methods)

    data(SpCas9, package="crisprBase", envir=environment())
    isCas9 <- .identicalNucleases(crisprNuclease, SpCas9)
    if (!is.null(fastaFile) & !is.null(chromatinFiles) & modality %in% c("CRISPRi", "CRISPRi") & isCas9){
        if (verbose){
            cat("[precomputeGuides] Adding CRISPRai scores \n")
        }
        out <- addCrispraiScores(out,
                                 gr=gr,
                                 tssObject=tssObject,
                                 modality=modality,
                                 chromatinFiles=chromatinFiles,
                                 fastaFile=fastaFile)
    }

    if (modality=='CRISPRko' & isCas9){
        if (verbose){
            cat("[precomputeGuides] Adding CFD scores annotation \n")
            out <- addOffTargetScores(out, max_mm=max_mm)
        }
    }
 
    if (!is.null(vcf)){
        if (verbose){
            cat("[precomputeGuides] Adding SNP annotation \n")
        }
        out <- addSNPAnnotation(out, vcf=vcf)
    }
    return(out)
}


