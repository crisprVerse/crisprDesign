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
#' @param commonExon Should gRNAs targeting common exons by prioritized?
#'     FALSE by default. If TRUE, \code{tx_id} must be provided.
#' @param modality String specifying the CRISPR modality. Should be one
#'     of the following: "CRISPRko", "CRISPRa", or "CRISPRi". 
#' @param useDistanceToTss Should distance to TSS be used to rank gRNAs 
#'     for CRISPRa and CRISPRi applications? TRUE by default. 
#'     For SpCas9 and human targets, this should be set to FALSE if
#'     \code{addCrispraiScores} was used. 
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
                        commonExon=FALSE,
                        modality=c("CRISPRko",
                                   "CRISPRa",
                                   "CRISPRi"),
                        useDistanceToTss=TRUE
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
                                              modality=modality)
    guideSet <- .getDefaultCompositeScores(guideSet)
    if (modality=="CRISPRko"){
        rankingCols <- c("round",
                         "score_exon",
                         "score_cds",
                         "score_conservation_binary",
                         "score_composite")
    } else if (modality=="CRISPRa" | modality=="CRISPRi"){
        if (useDistanceToTss){
            rankingCols <- c("round",
                             "distBin",
                             "score_composite")
        } else {
            rankingCols <- c("round",
                             "score_composite")
        }
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
                                          modality=c("CRISPRko",
                                                     "CRISPRa",
                                                     "CRISPRi")
){
    modality <- match.arg(modality)
    isKO <- modality=="CRISPRko"
    isA  <- modality=="CRISPRa"
    isI  <- modality=="CRISPRi"
    isAorI <- isA|isI
    
    cols <- colnames(mcols(guideSet))
    if (!"dist_to_tss" %in% cols  & isAorI){
        stop("dist_to_tss column is missing. Please call addDistanceToTss first.")
    }
    if (!"polyT" %in% cols){
        stop("poly column is missing. Please call addSequenceFeatures first.")
    }
    if (!"percentGC" %in% cols){
        stop("percentGC column is missing. Please call addSequenceFeatures first.")
    }
    n0     <- guideSet$n0
    n1     <- guideSet$n1
    n2     <- guideSet$n2
    n1_c   <- guideSet$n1_c
    n2_c   <- guideSet$n2_c
    n1_p   <- guideSet$n1_p
    n2_p   <- guideSet$n2_p
    polyT  <- guideSet$polyT
    hasSNP <- guideSet$hasSNP
    gc     <- guideSet$percentGC
  

    # Creating off-target bins:
    if (isAorI){
        guideSet$round <- NA
        guideSet$round[n0!=1 | (n1_p + n2_p)>5] <- 3
        guideSet$round[n0==1 & (n1_p + n2_p)<=5]  <- 2
        guideSet$round[n0==1 & n1_p==0 & n2_p==0] <- 1
    } else {
        guideSet$round <- NA
        guideSet$round[n0!=1 | (n1_c + n2_c)>5] <- 3
        guideSet$round[n0==1 & (n1_c + n2_c)<=5]  <- 2
        guideSet$round[n0==1 & n1_c==0 & n2_c==0] <- 1
    }
    guideSet$round[polyT] <- 4
    guideSet$round[gc<20] <- 4
    guideSet$round[gc>80] <- 4
    guideSet$round[is.na(n1)] <- 4
    guideSet$round[is.na(n2)] <- 4
    if (!is.null(hasSNP)){
        guideSet$round[hasSNP] <- 4
    }
 
    # Creating distance bins:
    if (isA){
        dist    <- guideSet$dist_to_tss
        distBin <- dist
        starts <- c(-150, -175, -200, -250, -300, -350, -400, -500)
        ends   <- c(-75,   -50,  -25,    0,    0,    0,    0,    0)
        nbins <- length(starts)
        for (k in rev(seq_len(nbins))){
            distBin[dist>=starts[k] & dist<=ends[k]] <- k
        }
        guideSet$distBin <- distBin
    } else if (isI){
        #bin1: 25 to 75
        #bin2: 0 to 25
        #bin3: 75 to 100 and 175 to 250
        #bin4: -25 to 0, and 100 to 175 and <250
        dist    <- guideSet$dist_to_tss
        distBin <- dist
        distBin[dist>=25 & dist<=75] <- 1
        distBin[dist>=0 & dist<25] <- 2
        distBin[(dist>75 & dist<=100) | (dist>=175 & dist<250)] <- 3
        distBin[dist<0 | dist >=250 | (dist>100 & dist<175)] <- 4
        guideSet$distBin <- distBin
    } else {
        guideSet$distBin <- NA
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
            scores <- "crispra"
        }
        if (hasCrispriScore){
            scores <- "crispri"
        }
        if (hasCrispraScore & hasCrispriScore){
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

