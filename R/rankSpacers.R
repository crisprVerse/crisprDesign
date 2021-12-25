#' @title Filter and rank spacers
#' @description Functions for filtering or ranking spacers within a
#'     \linkS4class{GuideSet} object according to supplied criteria.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param criteria A named list of values to filter or rank spacers
#'     List names must match names in \code{mcols(guideSet)} or
#'     column names of annotation list-columns, such as \code{geneAnnotation}.
#'     Duplicate names are permitted. See Details.
#' @param txId String specifying transcript ID. If \code{criteria} contains
#'     gene-level criteria, either \code{txId} or \code{geneId} must 
#'     be provided. See Details for gene-level criteria.
#' @param geneId String specifying gene ID. If \code{criteria} contains
#'     gene-level criteria, either \code{txId} or \code{geneId} must 
#'     be provided. See Details for gene-level criteria.
#' @param isoformAggFun String specifying the function name to be used 
#'     to aggregate gene-level data for gene-level \code{criteria} parameters
#'     when \code{is.null(txId)} and \code{!is.null(geneId)}.
#'     See Details for gene-level criteria.
#' 
#' @return \code{filterSpacers} - a \linkS4class{GuideSet} object filtered 
#'     by the values in \code{criteria}.
#' 
#' @return \code{rankSpacers} - a \linkS4class{GuideSet} object sorted by bins
#'     defined in \code{criteria}. A \code{rankings} list-column will also be
#'     stored in \code{mcols(GuideSet)}, which includes guide bin values for
#'     each element in \code{criteria}, and a \code{rank} column aggregated
#'     from all possible bin value combinations.
#'     
#' @details
#'     Use \code{validCriteria(guideSet)} for information on expected values
#'     for \code{criteria}, \code{txId}, and \code{geneId} arguments. There
#'     are four parameters:
#'     
#'     \itemize{
#'         \item \strong{attribute} — acceptable names for \code{criteria},
#'         which correspond to the column names in \code{mcols} or list-columns
#'         of \code{guideSet}.
#'         \item \strong{valueType} — type of values used to divide data for
#'         the specified attribute. There are four types, with specific
#'         requirements:
#'         \itemize{
#'             \item \bold{logical}: TRUE or FALSE. Indicates a preference for
#'             which value to retain (\code{filterSpacers}) or assign a higher
#'             rank (\code{rankSpacers}).
#'             \item \bold{asc}: a numeric vector in strictly ascending order
#'             (no duplicates). Sets binning intervals that are closed on the
#'             right.
#'             \item \bold{desc}: a numeric vector in strictly descending order
#'             (no duplicates). Sets binning intervals that are closed on the
#'             left.
#'             \item \bold{ranged}: a numeric vector of length 2 with values
#'             in the range 0-100 (duplicates permitted). Sets an inclusive
#'             percentile range for spacers to retain (\code{filterSpacers})
#'             or assign a higher rank (\code{rankSpacers}).
#'         }
#'         \item \strong{takesTxId} — Whether the attribute uses \code{txId}
#'         to obtain values for filtering or ranking.
#'         \item \strong{takesGeneId} — Whether the attribute uses
#'         \code{geneId} to obtain values for filtering or ranking.
#'     }
#'     
#'     If \code{criteria} contains an attribute for which only one of
#'     \code{takesTxId} or \code{takesGeneId} is \code{TRUE} then the
#'     corresponding must be provided. For some attributes, either ID is
#'     acceptable. In this case, the \code{txId} will be used if supplied,
#'     regardless of whether \code{geneId} is also supplied However, if
#'     \code{geneId} is supplied and \code{txId} is not, the values for that
#'     attribute will obtained via the \code{geneId} and aggregated according
#'     to \code{isoformAggFun}.
#'     
#'     \code{filterSpacers} filters for spacers meeting all supplied criteria,
#'     or, those having all attributes in \code{criteria} in the first (best)
#'     bin. Consequently, \code{asc} and \code{desc} values in \code{criteria}
#'     passed to \code{filterSpacers} need only define the limits of the first
#'     bin.
#'     
#'     \code{rankSpacers} assigns ranks to and sorts spacers according to
#'     \code{criteria}, with the order of attributes in \code{criteria}
#'     determining ranking priority. Attributes in \code{criteria} may be
#'     used multiple times.
#'     
#'      
#' @examples
#' \dontrun{
#' data(guideSetExample, package="crisprDesign")
#' guides <- guideSetExample
#' 
#' # annotate guides
#' guides <- addSequenceFeatures(guides)
#' guides <- addRestrictionEnzymes(guides)
#' guides <- addOnTargetScores(guides, methods=c('deephf', 'azimuth'))
#' guides <- addGeneAnnotation(guides, txObject=grListExample)
#' 
#' # check guide attributes available for filtering and ranking
#' validCriteria(guides)
#' 
#' filter_criteria <- list(
#'     EcoRI = FALSE,
#'     KpnI = FALSE,
#'     polyT = FALSE
#' )
#' guides <- filterSpacers(guides, filter_criteria)
#'   
#' rank_criteria <- list(
#'     isCommonCodingExon = TRUE,
#'     score_deephf = c(0.7, 0.6, 0.5),
#'     score_azimuth = c(0.7, 0.6, 0.5)
#' )
#' guides <- rankSpacers(guides, rank_criteria, geneId="ENSG00000120645")
#' guides
#' }
#' 
#' @author Luke Hoberecht
#' 
#' @name rankSpacers
NULL


#' @rdname rankSpacers
#' @export
filterSpacers <- function(guideSet,
                          criteria,
                          txId=NULL,
                          geneId=NULL,
                          isoformAggFun=c('max', 'mean', 'min')
){
    guideSet <- .validateGuideSet(guideSet) 
    .checkCriteriaIsNamedList(criteria)
    if (length(guideSet) == 0){
        return(guideSet)
    }
    isoformAggFun <- match.arg(isoformAggFun)
  
    guideRankings <- .rankSpacersOnCriteria(guideSet = guideSet,
                                            criteria = criteria,
                                            txId = txId,
                                            geneId = geneId,
                                            isoformAggFun = isoformAggFun)
  
    guideSet <- guideSet[guideRankings$rank == 1]
    return(guideSet)
}


#' @rdname rankSpacers
#' @export
rankSpacers <- function(guideSet,
                        criteria,
                        txId=NULL,
                        geneId=NULL,
                        isoformAggFun=c('max', 'mean', 'min')
){
    guideSet <- .validateGuideSet(guideSet) # why return value?
    .checkCriteriaIsNamedList(criteria)
    if (length(guideSet) == 0){
        return(guideSet)
    }
    isoformAggFun <- match.arg(isoformAggFun)
  
    guideRankings <- .rankSpacersOnCriteria(guideSet=guideSet,
                                            criteria=criteria,
                                            txId=txId,
                                            geneId=geneId,
                                            isoformAggFun=isoformAggFun)
  
    guideRankings <- guideRankings[order(guideRankings$rank), , drop=FALSE]
    guideSet <- guideSet[row.names(guideRankings)]
    guideSet$rankings <- guideRankings
    # apply ranking criteria to metadata(guideSet)?
    return(guideSet)
}


.checkCriteriaIsNamedList <- function(criteria){
    isList <- is.list(criteria)
    hasNames <- length(names(criteria)) > 0
    hasNoEmptyNames <- all(names(criteria) != '')
    hasNoNANames <- all(!is.na(names(criteria)))
    stopifnot('criteria must be a non-empty named list' = {
        isList && hasNames && hasNoEmptyNames && hasNoNANames
    })
}


.rankSpacersOnCriteria <- function(guideSet,
                                   criteria,
                                   txId,
                                   geneId,
                                   isoformAggFun
){
    criteriaParameters <- validCriteria(guideSet)
    binValues <- .binGuidesOnCriteria(guideSet,
                                      criteria,
                                      criteriaParameters=criteriaParameters,
                                      txId,
                                      geneId,
                                      isoformAggFun)
    binValuesAndRank <- .appendGuideRanks(
        criteriaParameters=criteriaParameters,
        criteria=criteria,
        binValues=binValues)
    return(binValuesAndRank)
}





#' @rdname rankSpacers
#' @importFrom utils stack
#' @export
validCriteria <- function(guideSet
){
    guideSet <- .validateGuideSet(guideSet) # why return value?
    validCriteria <- utils::stack(.validCriteriaByValueType(guideSet))
    colnames(validCriteria) <- c('attribute', 'valueType')
    validCriteria$valueType <- as.character(validCriteria$valueType)
    
    idUsage <- .validCriteriaByGenomicIdUse()
    validCriteria$takesTxId <- validCriteria$attribute %in% idUsage$txId
    validCriteria$takesGeneId <- validCriteria$attribute %in% idUsage$geneId
    
    guideSetColumns <- c(names(mcols(guideSet)),
                         names(geneAnnotation(guideSet)),
                         names(tssAnnotation(guideSet)),
                         names(enzymeAnnotation(guideSet))[-1])
    rowIndices <- validCriteria$attribute %in% guideSetColumns
    validCriteria <- validCriteria[rowIndices, , drop=FALSE]
    return(validCriteria)
}


#' @importFrom S4Vectors mcols
.validCriteriaByValueType <- function(guideSet
){
    list(
        logical=c('polyA', 'polyC', 'polyG', 'polyT', 'startingGGGGG',
                  'selfHairpin', 'backboneHairpin', 'cut_cds', 'cut_fiveUTRs',
                  'cut_threeUTRs', 'cut_introns', 'isCommonExon',
                  'isCommonCodingExon', 'hasSNP', 'inRepeats',
                  names(enzymeAnnotation(guideSet))[-1]), # exclude ID column
        asc=c('aminoAcidIndex', 'downstreamATG', 'percentCDS', 'percentTx',
              grep('^n[0-9](_[cp])?$', names(S4Vectors::mcols(guideSet)),
                   value=TRUE)),
        desc=c('nIsoforms', 'percentIsoforms', 'nCodingIsoforms',
               'percentCodingIsoforms',
               grep('^score_', names(S4Vectors::mcols(guideSet)),
                    value=TRUE)),
        ranged=c('percentGC')
    )
}


.validCriteriaByGenomicIdUse <- function(
){
    list(
        txId=c('cut_cds', 'cut_fiveUTRs', 'cut_threeUTRs', 'cut_introns',
               'aminoAcidIndex', 'downstreamATG', 'percentCDS', 'percentTx'),
        geneId=c('isCommonExon', 'isCommonCodingExon', 'aminoAcidIndex',
                 'downstreamATG', 'nIsoforms', 'percentIsoforms',
                 'nCodingIsoforms', 'percentCodingIsoforms', 'percentCDS',
                 'percentTx')
    )
}





#' @importFrom S4Vectors DataFrame
.binGuidesOnCriteria <- function(guideSet,
                                 criteria,
                                 criteriaParameters,
                                 txId,
                                 geneId,
                                 isoformAggFun
){
    rowIndices <- criteriaParameters$attribute %in% names(criteria)
    criteriaParameters <- criteriaParameters[rowIndices, , drop=FALSE]
    .checkBothGenomicIds(guideSet=guideSet,
                         criteriaParameters=criteriaParameters,
                         txId=txId,
                         geneId=geneId)
    binValues <- vapply(seq_along(criteria), function(x){
        criterionName <- names(criteria)[x]
        criterionValue <- criteria[[x]]
        rowIndex <- criteriaParameters$attribute == criterionName
        if (sum(rowIndex) == 0){
            stop('Criterion not found in guideSet: "', criterionName, '"',
                 .errorMessageHintAnnotation(criterionName))
        }
        criteriaParameters <- criteriaParameters[rowIndex, , drop=FALSE]
        .binGuidesOnCriterion(guideSet=guideSet,
                              criterionName=criterionName,
                              criterionValue=criterionValue,
                              criterionType=criteriaParameters$valueType,
                              takesTxId=criteriaParameters$takesTxId,
                              takesGeneId=criteriaParameters$takesGeneId,
                              txId=txId,
                              geneId=geneId,
                              isoformAggFun=isoformAggFun)
    }, FUN.VALUE = numeric(length(guideSet)))
  
    binValues <- S4Vectors::DataFrame(binValues,
                                      row.names=names(guideSet))
    colnames(binValues) <- .binValuesColnames(names(criteria))
    return(binValues)
}


.checkBothGenomicIds <- function(guideSet,
                                 criteriaParameters,
                                 txId,
                                 geneId
){
    .checkSingleGenomicId(guideSet=guideSet,
                          id='txId',
                          inputIdValue=txId,
                          criteriaUsingThisId=criteriaParameters$takesTxId,
                          criteriaUsingOtherId=criteriaParameters$takesGeneId)
    .checkSingleGenomicId(guideSet=guideSet,
                          id='geneId',
                          inputIdValue=geneId,
                          criteriaUsingThisId=criteriaParameters$takesGeneId,
                          criteriaUsingOtherId=criteriaParameters$takesTxId)
    takesEitherId <- any(criteriaParameters$takesTxId &
                             criteriaParameters$takesGeneId)
    noIdSupplied <- is.null(txId) && is.null(geneId)
    # only need to check if both are NULL, as both are separately checked above
    if (takesEitherId && noIdSupplied){
        stop('Supplied criteria requires either a txId or geneId.')
    }
    invisible(NULL)
}


.checkSingleGenomicId <- function(guideSet,
                                  id,
                                  inputIdValue,
                                  criteriaUsingThisId,
                                  criteriaUsingOtherId
){
    criteriaUsingThisIdOnly <- criteriaUsingThisId & !criteriaUsingOtherId
    criteriaUsingSuppliedId <- criteriaUsingThisId & !is.null(inputIdValue)
    needToCheckInputId <- any(criteriaUsingThisIdOnly |
                                  criteriaUsingSuppliedId)
    if (needToCheckInputId){
        if (is.null(inputIdValue)){
            stop('Supplied criteria requires a ', id)
        }
        isInvalidType <- !is.vector(inputIdValue, mode='character')
        hasInvalidLength <- length(inputIdValue) != 1
        if (isInvalidType || hasInvalidLength){
            stop(id, ' must be a length 1 character vector')
        }
        guideSetIdColname <- gsub('Id$', '_id', id)
        guideSetIds <- c(geneAnnotation(guideSet)[[guideSetIdColname]],
                         tssAnnotation(guideSet)[[guideSetIdColname]])
        guideSetIds <- unique(guideSetIds)
        guideSetIds <- guideSetIds[!is.na(guideSetIds)]
        if (!inputIdValue %in% guideSetIds){
            stop(id, ' not found in guideSet: ', inputIdValue)
        }
    }
    invisible(NULL)
}





.errorMessageHintAnnotation <- function(attribute
){
    annSources <- list(
        addGeneAnnotation=c(as.character(unlist(
            .validCriteriaByGenomicIdUse()))),
        addOffTargetScores=c('score_cfd', 'score_mit'),
        addOnTargetScores=c('score_azimuth', 'score_deephf', 'score_ruleset1',
                            'score_lindel', 'score_deepcpf1', 'score_enpamgb'),
        addPamScores=c('score_pam'),
        addRepeats=c('inRepeats'),
        addRestrictionEnzymes=c('EcoRI', 'KpnI', 'BsmBI', 'BsaI', 'BbsI',
                                'PacI', 'MluI'),
        addSequenceFeatures=c('percentGC', 'polyA', 'polyC', 'polyG', 'polyT',
                              'startingGGGGG', 'selfHairpin',
                              'backboneHairpin'),
        addSNPAnnotation=c('hasSNP'),
        addSpacerAlignments=c('n0', 'n0_c', 'n0_p', 'n1', 'n1_c', 'n1_p',
                              'n2', 'n2_c', 'n2_p', 'n3', 'n3_c', 'n3_p'),
        addTssAnnotation=c('dist_to_tss')
    )
    annSources <- stack(annSources)
    suggestedAnn <- annSources$ind[annSources$values == attribute]
    if (length(suggestedAnn) > 0){
        return(paste0('\nTry adding missing annotation with ',
                      suggestedAnn, '()\n'))
    }
    return(NULL)
}


.binGuidesOnCriterion <- function(guideSet,
                                  criterionName,
                                  criterionValue,
                                  criterionType,
                                  takesTxId,
                                  takesGeneId,
                                  txId,
                                  geneId,
                                  isoformAggFun
){
    valuesToBin <- .guideSetValuesToBin(guideSet=guideSet,
                                        criterionName=criterionName,
                                        takesTxId=takesTxId,
                                        takesGeneId=takesGeneId,
                                        txId=txId,
                                        geneId=geneId,
                                        isoformAggFun=isoformAggFun)
    binningFunction <- switch(criterionType,
                              logical=.binLogicalValues,
                              asc=.binAscValues,
                              desc=.binDescValues,
                              ranged=.binRangedValues
    )
    binningFunction(valuesToBin, criterionValue, criterionName)
}


#' @importFrom S4Vectors mcols
.guideSetValuesToBin <- function(guideSet,
                                 criterionName,
                                 takesTxId,
                                 takesGeneId,
                                 txId,
                                 geneId,
                                 isoformAggFun
){
    if (criterionName %in% names(S4Vectors::mcols(guideSet))){
        valuesToBin <- S4Vectors::mcols(guideSet)[[criterionName]]
    } else if (criterionName %in% names(enzymeAnnotation(guideSet))){
        valuesToBin <- enzymeAnnotation(guideSet)[[criterionName]]
    } else { # geneAnnotation
        valuesToBin <- .geneAnnotationValuesToBin(guideSet=guideSet,
                                                  attribute=criterionName,
                                                  takesTxId=takesTxId,
                                                  txId=txId,
                                                  takesGeneId=takesGeneId,
                                                  geneId=geneId,
                                                  isoformAggFun=isoformAggFun)
    }
    return(valuesToBin)
}


.geneAnnotationValuesToBin <- function(guideSet,
                                       attribute,
                                       takesTxId,
                                       txId,
                                       takesGeneId,
                                       geneId,
                                       isoformAggFun
){
    geneAnn <- geneAnnotation(guideSet, unlist=FALSE)
    getValuesUsingTxId <- takesTxId && !is.null(txId)
    valuesToBin <- lapply(geneAnn, function(x){
        if (getValuesUsingTxId){
            value <- x[[attribute]][x$tx_id == txId]
        } else {
            value <- x[[attribute]][x$gene_id == geneId]
            value <- value[!is.na(value)]
            if (takesTxId && takesGeneId && length(value) > 0){
                value <- switch(isoformAggFun,
                                max = max(value),
                                mean = mean(value),
                                min = min(value))
            } else {
                value <- unique(value)
            }
        }
        if (length(value) == 0){
            value <- NA
        }
        value
    })
    valuesToBin <- unlist(valuesToBin)
}


.binLogicalValues <- function(guideSetValues,
                              criterionValues,
                              criterionAttribute
){
    .checkLogicalCriterion(criterionValues, criterionAttribute)
    binValues <- factor(guideSetValues,
                        levels=c(criterionValues, !criterionValues))
    binValues <- as.numeric(binValues)
    binValues[is.na(binValues)] <- 3
    return(binValues)
}


.checkLogicalCriterion <- function(criterionValues,
                                   criterionAttribute
){
    hasWrongType <- !is.vector(criterionValues, mode='logical')
    hasBadLength <- length(criterionValues) != 1
    if (hasWrongType || hasBadLength){
        stop(criterionAttribute, ' must be a logical vector of length 1')
    }
    invisible(NULL)
}


.binAscValues <- function(guideSetValues,
                          criterionValues,
                          criterionAttribute
){
    .binNumericSeqValues(guideSetValues, criterionValues, criterionAttribute,
                         decreasing=FALSE)
}


.binDescValues <- function(guideSetValues,
                           criterionValues,
                           criterionAttribute
){
    .binNumericSeqValues(guideSetValues, criterionValues, criterionAttribute,
                         decreasing=TRUE)
}


.binNumericSeqValues <- function(guideSetValues,
                                 criterionValues,
                                 criterionAttribute,
                                 decreasing
){
    .checkNumericSeqCriterion(criterionValues, criterionAttribute, decreasing)
    criterionValues <- unique(c(criterionValues, Inf, -Inf))
    criterionValues <- sort(criterionValues, decreasing=decreasing)
    naReplacementValue <- criterionValues[length(criterionValues)]
    guideSetValues[is.na(guideSetValues)] <- naReplacementValue
    binValues <- cut(guideSetValues,
                     breaks=criterionValues,
                     labels=FALSE,
                     include.lowest=TRUE,
                     right=!decreasing)
    if (decreasing){
        binValues <- length(criterionValues) - binValues
    }
    return(binValues)
}


.checkNumericSeqCriterion <- function(criterionValues,
                                      criterionAttribute,
                                      decreasing
){
    hasWrongType <- !is.vector(criterionValues, mode='numeric')
    if (hasWrongType){
        stop(criterionAttribute, ' must be a numeric vector')
    }
    isNotStrictlyOrdered <- !identical(
        criterionValues,
        sort(unique(criterionValues), decreasing=decreasing)
    )
    if (isNotStrictlyOrdered){
        order <- ifelse(decreasing, 'descending', 'ascending')
        stop(criterionAttribute, ' must be a strictly ', order,
             ' numeric sequence')
    }
    invisible(NULL)
}


.binRangedValues <- function(guideSetValues,
                             criterionValues,
                             criterionAttribute
){
    .checkRangedCriterion(criterionValues, criterionAttribute)
    binValues <- guideSetValues >= criterionValues[1] &
        guideSetValues <= criterionValues[2]
    binValues <- factor(binValues, levels=c(TRUE, FALSE))
    binValues <- as.numeric(binValues)
    return(binValues)
}


.checkRangedCriterion <- function(criterionValues,
                                  criterionAttribute
){
    hasWrongType <- !is.vector(criterionValues, mode='numeric')
    hasBadLength <- length(criterionValues) != 2
    if (hasWrongType || hasBadLength){
        stop(criterionAttribute, ' must be a numeric vector of length 2')
    }
    hasBadRange <- min(criterionValues) < 0 || max(criterionValues) > 100
    if (hasBadRange){
        stop(criterionAttribute, ' must have values between 0-100')
    }
    isNotAscending <- !identical(criterionValues, sort(criterionValues))
    if (isNotAscending){
        stop(criterionAttribute, ' must be in ascending order')
    }
    invisible(NULL)
}





.binValuesColnames <- function(criteriaNames
){
    newColnames <- paste0('bin1_', criteriaNames)
    dup <- duplicated(newColnames)
    while(sum(dup) > 0){
        newColnames[dup] <- vapply(newColnames[dup], function(x){
            pattern <- '^bin([0-9]+)_'
            m <- regexec(pattern, x)
            currentBinNumber <- regmatches(x, m)[[1]][2]
            newBinNumber <- as.numeric(currentBinNumber) + 1
            newPrefix <- paste0('bin', newBinNumber, '_')
            gsub(pattern, newPrefix, x)
        }, FUN.VALUE = character(1))
        dup <- duplicated(newColnames)
    }
    return(newColnames)
}





.appendGuideRanks <- function(criteriaParameters,
                              criteria,
                              binValues
){
    binsToRankKey <- .binsToRankKey(criteriaParameters=criteriaParameters,
                                    criteria=criteria,
                                    binColumnNames=colnames(binValues))
    binValues$id <- row.names(binValues)
    mergedData <- merge(binValues, binsToRankKey)
    binValuesOrder <- match(binValues$id, mergedData$id)
    mergedData <- mergedData[binValuesOrder, , drop=FALSE]
    row.names(mergedData) <- mergedData$id
    newColOrder <- c('id', 'rank')
    newColOrder <- c(newColOrder, setdiff(colnames(mergedData), newColOrder))
    mergedData <- mergedData[, newColOrder, drop=FALSE]
    return(mergedData)
}


.binsToRankKey <- function(criteriaParameters,
                           criteria,
                           binColumnNames
){
    possibleBinValues <- lapply(seq_along(criteria), function(x){
        rowIndex <- criteriaParameters$attribute == names(criteria)[x]
        criterion <- criteriaParameters[rowIndex, , drop=FALSE]
        takesId <- criterion$takesTxId || criterion$takesGeneId
        seqLength <- switch(criterion$valueType,
            logical = ifelse(takesId, 3, 2),
            asc = length(unique(c(criteria[[x]], Inf, -Inf))) - 1,
            desc = length(unique(c(criteria[[x]], Inf, -Inf))) - 1,
            ranged = 2
        )
        seq_len(seqLength)
    })
    possibleBinValues <- rev(possibleBinValues)
    binsToRankKey <- expand.grid(possibleBinValues)
    binsToRankKey <- binsToRankKey[, rev(colnames(binsToRankKey)), drop=FALSE]
    colnames(binsToRankKey) <- binColumnNames
    binsToRankKey$rank <- seq_len(nrow(binsToRankKey))
    return(binsToRankKey)
}
