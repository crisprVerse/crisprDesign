# library(crisprDesign)
# library(crisprDesignGne)
# data(guideSetExampleFullAnnotation)
# guideSet <- guideSetExampleFullAnnotation
# tx_id = "ENST00000538872"
# conservationFile <- getConservationFiles()
# guideSet <- addConservationScores(guideSet,
#                                   conservationFile=conservationFile)
# guideSet <- rankSpacers(guideSet, 
#                         tx_id=tx_id)


#' @title Recommended gRNA ranking
#' @description Function for ranking spacers using 
#'     recommended crisprDesign criteria. CRISPRko, CRISPRa and CRISPRi
#'     modalities are supported. 
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param tx_id Optional string specifying transcript ID to use
#'     isoform-specific information for gRNA ranking.
#' @param useCodingInfo Should protein coding annotation be used to rank 
#'    gRNAs based on the number of off-targets? TRUE by default. 
#' @param commonExon Should gRNAs targeting common exons by prioritized?
#'     FALSE by default. If TRUE, \code{tx_id} must be provided.
#' @param modality String specifying the CRISPR modality. Should be one
#'     of the following: "CRISPRko", "CRISPRa", or "CRISPRi". 
#' 
#' @return A \linkS4class{GuideSet} object ranked from best to worst gRNAs,
#'     with a column \code{rank} stored in \code{mcols(guideSet)} indicating
#'     gRNA rank. 
#' 
#' @details 
#' 
#'     For each nuclease, we rank gRNAs based on several rounds of priority.
#'     For SpCas9, gRNAs with unique target sequences and without 
#'     1-or 2-mismatch off-targets located in coding regions are placed into 
#'     the first round. Then, gRNAs with a small number of one- or 
#'     two-mismatch off-targets  (less than 5) are placed into the second 
#'     round. Remaining gRNAs are placed into the third round. Finally, 
#'     any gRNAs overlapping a common SNP (human only), containing a polyT 
#'     stretch, or with extreme GC content (below 20% of above 80%) 
#'     are placed into the fourth round.
#' 
#'     If \code{tx_id} is specified, within each round of selection, gRNAs 
#'     targeting the first 85 percent of the specific transcript are 
#'     prioritized first. If \code{tx_id} is specified, and \code{commonExon}
#'     is set to TRUE, gRNAs targeting common exons across isoforms are also
#'     prioritized. If a conservation score is available, gRNAs
#'     targeting conserved regions (phyloP conservation score greater than 0), 
#'     are also prioritized.
#' 
#'     Within each bin, gRNAs are ranked by a composite 
#'     on-target activity rank to prioritize active gRNAs. The composite
#'     on-target activity rank is calculated by taking the average rank
#'     across the DeepHF and DeepSpCas9 scores for CRISPRko.
#'     For CRISPRa or CRISPRi, the CRISPRai scores are used if available.
#' 
#'     The process is identical for enAsCas12a, with the exception that the 
#'     enPAMGb method is used as the composite score.
#' 
#'     For CasRx, gRNAs targeting all isoforms of a given gene, with no 1- or 
#'     2-mismatch off-targets, are placed into the first round. gRNAs
#'     targeting at least 50 percent of the isoforms of a given gene, with 
#'     no 1- or 2-mismatch off-targets, are placed into the second round. 
#'     Remaining gRNAs are placed into the third round. Within each
#'     round of selection, gRNAs are further ranked by the CasRxRF
#'     on-target score. 
#' 
#' @examples
#' data(guideSetExampleFullAnnotation, package="crisprDesign")
#' gs <- rankSpacers(guideSetExampleFullAnnotation,
#'                   tx_id = "ENST00000538872")
#' gs
#' 
#' @author Luke Hoberecht, Jean-Philippe Fortin
#' 
#' @export
rankSpacers <- function(guideSet,
                        tx_id=NULL,
                        useCodingInfo=TRUE,
                        commonExon=FALSE,
                        modality=c("CRISPRko", "CRISPRa", "CRISPRi")
){
    modality <- match.arg(modality)
    crisprNuclease <- crisprNuclease(guideSet)
    data(SpCas9,     package="crisprBase", envir=environment())
    data(enAsCas12a, package="crisprBase", envir=environment())
    data(CasRx,      package="crisprBase", envir=environment())
    isCas9   <- .identicalNucleases(crisprNuclease, SpCas9)
    isCas12a <- .identicalNucleases(crisprNuclease, enAsCas12a)
    isCas13d <- .identicalNucleases(crisprNuclease, CasRx)

    # Adding isoform specific information:
    if (!is.null(tx_id)){
        guideSet <- addIsoformAnnotation(guideSet,
                                         tx_id=tx_id)
    }

    metacols <- mcols(guideSet)
    cols <- colnames(metacols)


    if ("percentCDS" %in% cols){
        guideSet$percentCDS[is.na(guideSet$percentCDS)] <- 100
        guideSet$score_cds <- ifelse(guideSet$percentCDS<=85,1,0)
    }
    if ("score_conservation" %in% cols){
        newScore <- ifelse(guideSet$score_conservation>=0,1,0)
        guideSet$score_conservation_binary <- newScore
    }
    if ("isCommonCodingExon" %in% cols){
        exon <- guideSet$isCommonCodingExon
        exon <- ifelse(exon,1,0)
        guideSet$score_exon <- exon
    }
    guideSet <- .createDefaultSelectionRounds(guideSet,
                                              useCodingInfo=useCodingInfo)
    guideSet <- .getDefaultCompositeScores(guideSet)
    if (modality=="CRISPRko"){
        rankingCols <- c("round",
                         "score_exon",
                         "score_cds",
                         "score_conservation_binary",
                         "score_composite")
    } else if (modality=="CRISPRa" | modality=="CRISPRi"){
        rankingCols <- c("round",
                         "score_composite")
    }

    rankingCols <- intersect(rankingCols, colnames(mcols(guideSet)))
    ranks <- lapply(rankingCols, function(x) mcols(guideSet)[[x]])
    names(ranks) <- rankingCols

    # Transforming scores for ascending ranking:
    scoreCols <- which(grepl("score", rankingCols))
    if (length(scoreCols)>0){
        ranks[scoreCols] <- lapply(ranks[scoreCols], function(col){
            -col
        })
    }
    guideSet <- guideSet[do.call("order", ranks)]
    guideSet$rank <- seq_len(length(guideSet))
    return(guideSet)
}



.createDefaultSelectionRounds <- function(guideSet,
                                          useCodingInfo=TRUE
){
    n0     <- guideSet$n0
    n1     <- guideSet$n1
    n2     <- guideSet$n2
    n1_c   <- guideSet$n1_c
    n2_c   <- guideSet$n2_c
    polyT  <- guideSet$polyT
    hasSNP <- guideSet$hasSNP
    gc     <- guideSet$percentGC

    # Creating bins:
    guideSet$round <- NA
    if (useCodingInfo){
        if (is.null(n1_c) | is.null(n2_c)){
            stop("n1_c or n2_c columns don't exist. Please", 
                 " use addSpacerAlignments first with a specified",
                 " txObject.")
        }
        guideSet$round[n0!=1 | (n1 + n2_c)>5] <- 3
        guideSet$round[n0==1 & (n1_c + n2_c)<=5]  <- 2
        guideSet$round[n0==1 & n1_c==0 & n2_c==0] <- 1
    } else {
        guideSet$round[n0!=1 | n1>5]  <- 3
        guideSet$round[n0==1 & n1<=5] <- 2
        guideSet$round[n0==1 & n1==0 & n2<=3] <- 1
    }

    # Undesirable gRNAs:
    if (is.null(polyT) | is.null(gc)){
        stop("poly or percentGC columns are missing. Please use",
             " addSequenceFeatures first.")
    }
    guideSet$round[polyT] <- 4
    guideSet$round[gc<20] <- 4
    guideSet$round[gc>80] <- 4
    guideSet$round[is.na(n1)] <- 4
    guideSet$round[is.na(n2)] <- 4
    if (!is.null(hasSNP)){
        guideSet$round[hasSNP] <- 4
    }
    return(guideSet)
}


.getDefaultCompositeScores <- function(guideSet){
    crisprNuclease <- crisprNuclease(guideSet)
    data(SpCas9,     package="crisprBase", envir=environment())
    data(enAsCas12a, package="crisprBase", envir=environment())
    data(CasRx,      package="crisprBase", envir=environment())
    isCas9   <- .identicalNucleases(crisprNuclease, SpCas9)
    isCas12a <- .identicalNucleases(crisprNuclease, enAsCas12a)
    isCas13d <- .identicalNucleases(crisprNuclease, CasRx)
    if (isCas9){
        scores <- c("deephf", "deepspcas9")
        hasCrispraScore <- "score_crispra" %in% colnames(mcols(guideSet))
        hasCrispriScore <- "score_crispri" %in% colnames(mcols(guideSet))
        if (hasCrispraScore){
            scores <- "score_crispra"
        }
        if (hasCrispriScore){
            scores <- "score_crispri"
        }
        if (hasCrispriScore & hasCrispriScore){
            stop("Both CRISPRi and CRISPRa scores are added.",
                 "Please remove one.")
        }
    } else if (isCas12a){
        scores <- c("enpamgb")
    } else if (isCas13d){
        scores <- c("casrxrf")
    }
    scoreCols <- paste0("score_",scores)

    cols <- colnames(mcols(guideSet))
    missing <- setdiff(scoreCols, cols)
    if (length(missing)>0){
        missing <- paste0(missing, collapse=",")
        stop("The following on-target scoring methods are missing ",
             " and are needed for default rankings: ", missing)   
    }
    guideSet <- addCompositeScores(guideSet, methods=scores)
    return(guideSet)
}

