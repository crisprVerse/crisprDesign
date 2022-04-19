#' @title Add CRISPRa/CRISPRi on-target scores to a \linkS4class{GuideSet} object.
#' @description Add CRISPRa/CRISPRi on-target scores to a
#'    \linkS4class{GuideSet} object. Only available for SpCas9, and for 
#'    hg38 genome. Requires \pkg{crisprScore} package to be installed.
#' 
#' @param guideSet A \linkS4class{GuideSet} object. 
#' @param gr A \linkS4class{GRanges} object derived from \code{queryTss} used
#'     to produce the \code{guideSet} object.
#' @param tssObject  A \linkS4class{GRanges} object containing TSS coordinates
#'     and annotation. The following columns must be present:
#'     "ID", promoter", "tx_id" and "gene_symbol".
#'     
#' @param modality String specifying which modality is used.
#'     Must be either "CRISPRi" or "CRISPRa". 
#' @param chromatinFiles Named character vector of length 3 specifying
#'     BigWig files containing chromatin accessibility data. See
#'     crisprScore vignette for more information.
#' @param fastaFile String specifying fasta file of the hg38 genome.
#' 
#' @return \code{guideSet} with an added column for the CRISPRai score.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addOnTargetScores}} to add other on-target scores.
#' 
#' @export
#' @importFrom crisprScore getCrispraiScores
addCrispraiScores <- function(guideSet,
                              gr,
                              tssObject,
                              modality=c("CRISPRi", "CRISPRa"),
                              chromatinFiles=NULL,
                              fastaFile=NULL
){
    crisprNuclease <- crisprNuclease(guideSet)
    data(SpCas9,
         package="crisprBase",
         envir=environment())
    if (!.identicalNucleases(crisprNuclease, SpCas9)){
        stop("[addCrispraiScores] Only SpCas9 is supported at the moment.")
    }

    if (genome(guideSet)[1]!="hg38"){
        stop("addCrispraiScores] Only hg38 genome supported at the moment.")        
    }
    modality <- match.arg(modality)
    tssFrame  <- .prepareTssFrame(tssObject)
    grnaFrame <- .prepareGrnaFrame(guideSet, gr)
    scores <- crisprScore::getCrispraiScores(sgrna_df=grnaFrame,
                                             tss_df=tssFrame,
                                             chromatinFiles=chromatinFiles,
                                             fastaFile=fastaFile,
                                             modality=modality)
    scores <- scores[match(names(guideSet), rownames(scores)),1]
    if (modality=="CRISPRa"){
        mcols(guideSet)$score_crispra <- scores
    } else {
        mcols(guideSet)$score_crispri <- scores
    }
    return(guideSet)
}



.prepareTssFrame <- function(tssObject){
    tssObject <- as.data.frame(tssObject)

    cols <- c("ID", "gene_symbol", "promoter","tx_id")
    if (!all(cols %in% colnames(tssObject))){
        choices <- setdiff(cols, colnames(tssObject))
        stop("The following columns are missing in the tssObject: \n \t",
             paste0(choices, collapse=", "),".")
    }

    out <- data.frame(tss_id=tssObject$ID,
                      gene_symbol=tssObject$gene_symbol,
                      promoter=tssObject$promoter,
                      transcripts=tssObject$tx_id,
                      position=tssObject$start,
                      strand=tssObject$strand,
                      chr=tssObject$seqnames)

    # Check if there are any missing values:
    if (sum(is.na(out$gene_symbol))>0){
        stop("gene_symbol has some missing values.")
    }
    if (sum(is.na(out$promoter))>0){
        stop("promoter has some missing values.")
    }
    #if (sum(is.na(out$transcripts))>0){
    #    stop("tx_id has some missing values.")
    #}
    if (sum(is.na(out$position))>0){
        stop("start has some missing values.")
    }
    if (sum(is.na(out$strand))>0){
        stop("strand has some missing values.")
    }
    if (sum(is.na(out$seqnames))>0){
        stop("strand has some missing values.")
    }
    return(out)
}




.prepareGrnaFrame <-function(guideSet, gr){

    len <- spacerLength(guideSet)
    seqs <- spacers(guideSet, as.character=TRUE)
    if (len!=19){
        if (len==20){
            seqs <- substr(seqs, 2,20)
        } else{
            stop("spacer length must be of length 19 or 20.")
        }
    } 

    cols <- c("ID")
    colnames <- colnames(mcols(gr))
    if (!all(cols %in% colnames)){
        choices <- setdiff(cols, colnames)
        stop("The following columns are missing in the gr object: \n \t",
             paste0(choices, collapse=", "),".")
    }


    ids <- gr$ID[match(guideSet$region, names(gr))]
    if (sum(is.na(ids))>0){
        stop("Some of the guideSet regions cannot be found in the gr object.")
    }
    out <- data.frame(grna_id=names(guideSet),
                      tss_id=ids, 
                      pam_site=pamSites(guideSet),
                      strand=as.character(strand(guideSet)),
                      spacer_19mer=seqs)
    return(out)
}







# library(crisprDesign)
# library(crisprDesignGne)
# bsgenome <- getGenomePackage()
# gr <- queryTss(tssObjectExample,
#                tss_window=c(-1000,1000),
#                queryColumn="gene_symbol",
#                queryValue="IQSEC3")
# guideSet <- findSpacers(gr, bsgenome=bsgenome)
# chromatinFiles <- getChromatinFiles()
# fastaFile <- getGenomeFasta()


# tssObject=tssObjectExample


# scores_a <- addCrispraiScores(guideSet,
#                               gr,
#                               tssObject,
#                               "CRISPRa",
#                               chromatinFiles=chromatinFiles,
#                               fastaFile=fastaFile)
# scores_i <- addCrispraiScores(guideSet,
#                               gr,
#                               tssObject,
#                               "CRISPRi",
#                               chromatinFiles=chromatinFiles,
#                               fastaFile=fastaFile)

# # For CRISPRa:
# par(mfrow=c(1,2))
# wh <- which(guideSet$region=="region_1")
# plot(pamSites(guideSet)[wh]-66767, scores_a[wh])
# abline(v=-150)
# wh <- which(guideSet$region=="region_2")
# plot(pamSites(guideSet)[wh]-77376, scores_a[wh])
# abline(v=-150)


# # For CRISPRi:
# par(mfrow=c(1,2))
# wh <- which(guideSet$region=="region_1")
# plot(pamSites(guideSet)[wh]-66767, scores_i[wh])
# abline(v=0)
# wh <- which(guideSet$region=="region_2")
# plot(pamSites(guideSet)[wh]-77376, scores_i[wh])
# abline(v=0)

















