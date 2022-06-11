#' @title Add on-target scores to a \linkS4class{GuideSet} object.
#' @description Add on-target scores to a \linkS4class{GuideSet} object
#'    for all methods available in the \pkg{crisprScore} package for a given
#'    CRISPR nuclease. Requires \pkg{crisprScore} package to be installed.
#' 
#' @param guideSet A \linkS4class{GuideSet} object. 
#' @param enzyme Character string specifying the Cas9 variant to be used
#'     for DeepHF scoring. Wildtype Cas9 (WT) by default. See details below.
#' @param promoter Character string speciyfing promoter used for expressing 
#'     sgRNAs for wildtype Cas9 (must be either "U6" or "T7") for DeepHF
#'     scoring. "U6" by default. 
#' @param methods Character vector specifying method names for on-target
#'     efficiency prediction algorithms.
#' 
#' @return \code{guideSet} with columns of on-target scores appended in
#'     \code{mcols(guideSet)}.
#' 
#' @details See \pkg{crisprScore} package for a description of each score.
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{\link{addOffTargetScores}} to add off-target scores.
#' 
#' @examples
#' if (interactive()){
#'     gs <- findSpacers("CCAACATAGTGAAACCACGTCTCTATAAAGAATAAAAAATTAGCCGGGTTA")
#'     gs <- addOnTargetScores(gs)
#' }
#' 
#' @export
addOnTargetScores <- function(guideSet,
                              enzyme=c("WT", "ESP", "HF"),
                              promoter=c("U6", "T7"),
                              methods=c("azimuth",
                                        "deephf",
                                        "ruleset1",
                                        "lindel",
                                        "deepcpf1",
                                        "enpamgb",
                                        "crisprscan",
                                        "casrxrf")
){
    guideSet <- .validateGuideSetOrPairedGuideSet(guideSet)
    if (.isGuideSet(guideSet)){
        out <- addOnTargetScores_guideset(guideSet,
                                          enzyme=enzyme,
                                          promoter=promoter,
                                          methods=methods)
    } else if (.isPairedGuideSet(guideSet)){
        unifiedGuideSet <- .pairedGuideSet2GuideSet(guideSet)
        unifiedGuideSet <- addOnTargetScores_guideset(unifiedGuideSet,
                                                      enzyme=enzyme,
                                                      promoter=promoter,
                                                      methods=methods)
        out <- .addColumnsFromUnifiedGuideSet(guideSet,
                                              unifiedGuideSet)
    }
    return(out)
}




#' @importFrom S4Vectors mcols<-
addOnTargetScores_guideset <- function(guideSet,
                                       enzyme=c("WT", "ESP", "HF"),
                                       promoter=c("U6", "T7"),
                                       methods=c("azimuth",
                                                 "deephf",
                                                 "ruleset1",
                                                 "lindel",
                                                 "deepcpf1",
                                                 "enpamgb",
                                                 "crisprscan",
                                                 "casrxrf")
){
    guideSet <- .validateGuideSet(guideSet)
    enzyme <- match.arg(enzyme)
    promoter <- match.arg(promoter)
    crisprNuclease <- crisprNuclease(guideSet)
    
    methods <- .validateOnTargetScoreMethods(methods=methods,
                                             crisprNuclease=crisprNuclease)

    valid <- .validSpacersForOnTargetScores(guideSet=guideSet,
                                            crisprNuclease=crisprNuclease)
    if (!any(valid)){
        warning("No valid chromosome data or PAM sequences; ",
                "cannot calculate on-target scores.",
                immediate.=TRUE)
    }
    for (i in methods){
        if (any(valid)){
            status <- paste0("[addOnTargetScores] Adding ", i,
                             " scores. \n")
            message(status)
        }
        scores <- .getOnTargetScores(guideSet=guideSet[valid],
                                     method=i,
                                     promoter=promoter,
                                     enzyme=enzyme)
        scoreColname <- paste0("score_", i)
        S4Vectors::mcols(guideSet)[[scoreColname]] <- rep(NA,
                                                          length(guideSet))
        S4Vectors::mcols(guideSet)[[scoreColname]][valid] <- scores
    }
    return(guideSet)
}






#' @importFrom utils data
.validateOnTargetScoreMethods <- function(methods,
                                          crisprNuclease
){
    utils::data("scoringMethodsInfo",
                package="crisprScore",
                envir=environment())
    data(SpCas9, package="crisprBase", envir=environment())
    data(AsCas12a, package="crisprBase", envir=environment())
    data(enAsCas12a, package="crisprBase", envir=environment())
    data(CasRx, package="crisprBase", envir=environment())
    if (.identicalNucleases(crisprNuclease, SpCas9)){
        choices <- c("azimuth",
                     "deephf",
                     "lindel",
                     "ruleset1",
                     "crisprscan")
    } else if (.identicalNucleases(crisprNuclease, AsCas12a)){
        choices <- c("deepcpf1")
    } else if (.identicalNucleases(crisprNuclease, enAsCas12a)){
        choices <- c("deepcpf1", "enpamgb")
    } else if (.identicalNucleases(crisprNuclease, CasRx)){
        choices <- c("casrxrf")
    } else {
        stop("No scoring method found for crisprNuclease \n")
    }
    badMethods <- setdiff(methods, scoringMethodsInfo$method)
    #badMethods <- setdiff(badMethods, "casrxrf")
    if (length(badMethods) > 0){
        stop("Scoring methods not recognized: ",
             paste(badMethods, collapse=", "))
    }
    ## add message for specified methods to be ignored?
    methods <- intersect(choices, methods)
    if (length(methods) == 0){
        stop("No valid methods specified for the crisprNuclease\n")
    }
    return(methods)
}



#' @importFrom stats complete.cases
.validSpacersForOnTargetScores <- function(guideSet,
                                           crisprNuclease
){
    df <- as.data.frame(guideSet)
    cols  <- c("seqnames", "start", "end", "pam")
    valid <- stats::complete.cases(df[, cols, drop=FALSE]) 
    # For Cas9, scores only work with canonical NGG PAM site:
    data(SpCas9, package="crisprBase", envir=environment())
    if (isTRUE(all.equal(crisprNuclease, SpCas9))){
        pams <- c("AGG", "TGG", "CGG", "GGG")
        validPam <- df[["pam"]] %in% pams
        valid <- valid & validPam
    }
    return(valid)
}





#' @author Jean-Philippe Fortin
#' 
#' @importFrom utils data
#' @importFrom crisprScore getDeepHFScores
#' @importFrom crisprScore getAzimuthScores
#' @importFrom crisprScore getRuleSet1Scores
#' @importFrom crisprScore getDeepCpf1Scores
#' @importFrom crisprScore getLindelScores
#' @importFrom crisprScore getEnPAMGBScores
#' @importFrom crisprScore getCRISPRscanScores
.getOnTargetScores <- function(guideSet,
                               method,
                               enzyme,
                               promoter
){
    if (method=="casrxrf"){
        scores <- .getCasRxRFScores(guideSet)
    } else {
        utils::data("scoringMethodsInfo",
                    package="crisprScore",
                    envir=environment())
        roster <- scoringMethodsInfo
        roster <- roster[roster$method == method, , drop=FALSE]
        left  <- roster$left
        right <- roster$right
        extendedSequences <- .getExtendedSequences(guideSet,
                                                   start=left,
                                                   end=right)
        good <- !is.na(extendedSequences)
        scores <- rep(NA, length(extendedSequences))
        if (any(good)){
            seqs <- extendedSequences[good]
            if (method == "deephf"){
                results <- crisprScore::getDeepHFScores(seqs,
                                                        enzyme=enzyme,
                                                        promoter=promoter)
            } else {
              scoreFun <- switch(method,
                                 "azimuth"=crisprScore::getAzimuthScores,
                                 "ruleset1"=crisprScore::getRuleSet1Scores,
                                 "deepcpf1"=crisprScore::getDeepCpf1Scores,
                                 "lindel"=crisprScore::getLindelScores,
                                 "enpamgb"=crisprScore::getEnPAMGBScores,
                                 "crisprscan"=crisprScore::getCRISPRscanScores)
              results <- scoreFun(seqs)
            }
            scores[good] <- results$score
        }
    } 
    return(scores)
}



#' @importFrom crisprScore getCasRxRFScores
.getCasRxRFScores <- function(guideSet){
    spacerLen <- spacerLength(guideSet)
    mrnaSequence <- metadata(guideSet)$customSequences
    if (length(mrnaSequence)>1){
        stop("mrnaSequence must be of length 1 for CasRxRF scoring")
    }
    if (spacerLen != 23){
        stop("Spacer length must be 23 to use CasRxRF")
    }
    scores <- crisprScore::getCasRxRFScores(mrnaSequence=mrnaSequence)
    wh <- match(spacers(guideSet, as.character=TRUE), scores$spacer)
    scores <- scores[wh,,drop=FALSE]
    out <- scores[["standardizedScore"]]
    return(out)
}








