#' @title Add on-target scores to a \linkS4class{GuideSet} object.
#' @description Add on-target scores to a \linkS4class{GuideSet} object
#'    for all methods available in the \pkg{crisprScore} package for a given
#'    CRISPR nuclease. Requires \pkg{crisprScore} package to be installed.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param enzyme Character string specifying the Cas9 variant to be used
#'     for DeepHF scoring. Wildtype Cas9 (WT) by default. See details below.
#' @param promoter Character string speciyfing promoter used for expressing 
#'     sgRNAs for wildtype Cas9 (must be either "U6" or "T7") for DeepHF
#'     scoring. "U6" by default. 
#' @param tracrRNA String specifying which tracrRNA is used for SpCas9
#'     Must be either "Hsu2013" (default) or "Chen2013". Only used for
#'     the RuleSet3 method.
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
#' @rdname addOnTargetScores
setMethod("addOnTargetScores", "GuideSet",
    function(object,
             enzyme=c("WT", "ESP", "HF"),
             promoter=c("U6", "T7"),
             tracrRNA=c("Hsu2013","Chen2013"),
             methods=c("azimuth",
                       "ruleset1",
                       "ruleset3",
                       "lindel",
                       "deepcpf1",
                       "deephf",
                       "deepspcas9",
                       "enpamgb",
                       "casrxrf",
                       "crisprater",
                       "crisprscan")
){
    object <- .validateGuideSet(object)
    enzyme <- match.arg(enzyme)
    promoter <- match.arg(promoter)
    tracrRNA <- match.arg(tracrRNA)
    crisprNuclease <- crisprNuclease(object)
    
    methods <- .validateOnTargetScoreMethods(methods=methods,
                                             crisprNuclease=crisprNuclease)

    valid <- .validSpacersForOnTargetScores(guideSet=object,
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
        scores <- .getOnTargetScores(guideSet=object[valid],
                                     method=i,
                                     promoter=promoter,
                                     tracrRNA=tracrRNA,
                                     enzyme=enzyme)
        scoreColname <- paste0("score_", i)
        S4Vectors::mcols(object)[[scoreColname]] <- rep(NA,
                                                          length(object))
        S4Vectors::mcols(object)[[scoreColname]][valid] <- scores
    }
    return(object)
})




#' @rdname addOnTargetScores
#' @export
setMethod("addOnTargetScores", "PairedGuideSet",
          function(object,
                   enzyme=c("WT", "ESP", "HF"),
                   promoter=c("U6", "T7"),
                   tracrRNA=c("Hsu2013","Chen2013"),
                   methods=c("azimuth",
                             "ruleset1",
                             "ruleset3",
                             "lindel",
                             "deepcpf1",
                             "deephf",
                             "deepspcas9",
                             "enpamgb",
                             "crisprater",
                             "crisprscan",
                             "casrxrf")
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addOnTargetScores(unifiedGuideSet,
                                         enzyme=enzyme,
                                         promoter=promoter,
                                         tracrRNA=tracrRNA,
                                         methods=methods)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addOnTargetScores
#' @export
setMethod("addOnTargetScores", "NULL", function(object){
    return(NULL)
})





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
                     "deepspcas9",
                     "lindel",
                     "ruleset1",
                     "ruleset3",
                     "crisprater",
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
#' @importFrom crisprScore getDeepSpCas9Scores 
#' @importFrom crisprScore getAzimuthScores
#' @importFrom crisprScore getRuleSet1Scores
#' @importFrom crisprScore getRuleSet3Scores
#' @importFrom crisprScore getDeepCpf1Scores
#' @importFrom crisprScore getLindelScores
#' @importFrom crisprScore getEnPAMGBScores
#' @importFrom crisprScore getCRISPRscanScores
#' @importFrom crisprScore getCRISPRaterScores
.getOnTargetScores <- function(guideSet,
                               method,
                               enzyme,
                               promoter,
                               tracrRNA
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
            } else if (method == "ruleset3"){
                results <- crisprScore::getRuleSet3Scores(seqs,
                                                          tracrRNA=tracrRNA)
            } else {
              scoreFun <- switch(method,
                                 "azimuth"=crisprScore::getAzimuthScores,
                                 "ruleset1"=crisprScore::getRuleSet1Scores,
                                 "deepcpf1"=crisprScore::getDeepCpf1Scores,
                                 "deepspcas9"=crisprScore::getDeepSpCas9Scores,
                                 "lindel"=crisprScore::getLindelScores,
                                 "enpamgb"=crisprScore::getEnPAMGBScores,
                                 "crisprater"=crisprScore::getCRISPRaterScores,
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








