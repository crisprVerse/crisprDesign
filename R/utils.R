#' @importFrom methods is
NULL

#' @importFrom utils globalVariables
utils::globalVariables(c("restrictionEnzymes",
                         "scoringMethodsInfo"))
utils::globalVariables(c("SpCas9",
                         "AsCas12a",
                         "enAsCas12a",
                         "CasRx"))


STOP_CODONS <- c("TAG","TAA","TGA")


#' @importFrom GenomeInfoDb seqnames
#' @export
GenomeInfoDb::seqnames


#' @importFrom S4Vectors mcols
#' @export
S4Vectors::mcols


.default_tss_window <- c(-500, 500)


.isGRanges <- function(x){
    methods::is(x, "GRanges")
}


.validateGRanges <- function(obj){
    if (!.isGRanges(obj)){
        stop("Object must be a GRanges.")
    } 
    return(obj)
}


.isGRangesList <- function(x){
    methods::is(x, "GRangesList")
}


.validateGRangesList <- function(obj){
    if (.isTxDb(obj)){
        return(TxDb2GRangesList(obj))
    }
    stopifnot("Object must be a TxDb or GRangesList." = {
        .isGRangesList(obj)
    })
    fields <- c("transcripts",
                "exons",
                "cds",
                "fiveUTRs",
                "threeUTRs",
                "introns",
                "tss")
    if (!all(fields %in% names(obj))){
        stop("GRangesList is missing required genomic regions. ",
             "Use TxDb2GRangesList to get properly formatted GRangesList")
    }
    return(obj)
}


.getTx2GeneTable <- function(txObject){
    cols <- c("gene_id", "tx_id", "gene_symbol")
    tx2Gene <- mcols(txObject[["exons"]])[, cols,drop=FALSE]
    tx2Gene <- as.data.frame(tx2Gene)
    tx2Gene <- tx2Gene[!duplicated(tx2Gene),]
    rownames(tx2Gene) <- NULL
    return(tx2Gene)
}


.isTxDb <- function(x){
    methods::is(x, "TxDb")
}


.isDNAStringSet <- function(x){
    methods::is(x, "DNAStringSet")
}


.isGuideSet <- function(object){
    methods::is(object, "GuideSet")
}


.validateGuideSet <- function(object){
    if (!.isGuideSet(object)){
        stop("Object must be a GuideSet")
    }
    return(object)
}


.validateGuideSetSequences <- function(argument,
                                       value){
    if (any(value == "")){
        stop(sprintf("'%s' cannot contain empty strings", argument))
    }
    return(value)
}


#' @importFrom methods is
.validateCrisprNuclease <- function(crisprNuclease){
    if (is.null(crisprNuclease)){
        crisprNuclease <- .getDefaultCrisprNuclease()
    }
    if (!methods::is(crisprNuclease, "CrisprNuclease")){
        stop("Provided nuclease must be a 'CrisprNuclease' object.")
    }
    return(crisprNuclease)
}


#' @importFrom utils data
.getDefaultCrisprNuclease <- function(type=c("Cas9", "Cas12a")){
    type <- match.arg(type)
    nuc <- switch(type,
                  "Cas9"="SpCas9",
                  "Cas12a"="AsCas12a")
    utils::data(list=nuc,
                package="crisprBase",
                envir=environment())
    nuc <- get(nuc)
    return(nuc)
}


.isBSGenome <- function(bsgenome){
    if (!methods::is(bsgenome, "BSgenome")){
        stop("Provided genome must be a 'BSgenome' object. ")
    }
}


#' @importFrom GenomeInfoDb genome
.getGenome <- function(guideSet){
    genome <- unique(GenomeInfoDb::genome(guideSet))
    stopifnot("Multiple genomes found for GuideSet object" = {
        length(genome) == 1
    })
    return(genome)
}


#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics rownames
.asDataFrame <- function(gr){
    spacerIds <- names(gr)
    # intermediate data.frame flattens gr
    df <- as.data.frame(gr)
    colnames(df)[colnames(df) == "seqnames"] <- "chr"
    colnames(df)[colnames(df) == "pos"] <- "anchor_site"
    df <- S4Vectors::DataFrame(df)
    BiocGenerics::rownames(df) <- spacerIds
    return(df)
}


#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics unlist
.isAlignmentsColumn <- function(guideSet,
                                columnName
){
    column <- S4Vectors::mcols(guideSet)[[columnName]]
    columnNameNotFound <- is.null(column)
    columnIsNotGRangesList <- !.isGRangesList(column)
    if (columnNameNotFound || columnIsNotGRangesList){
        return(FALSE)
    }
    column <- BiocGenerics::unlist(column)
    alignmentColumnNames <- c("spacer", "protospacer", "pam", "pam_site",
                              "n_mismatches", "canonical", "cut_site")
    hasAlignmentsColumns <- all(alignmentColumnNames %in%
                                    colnames(S4Vectors::mcols(column)))
    return(hasAlignmentsColumns)
}


#' @importFrom S4Vectors mcols
.validateAnchor <- function(anchor, gr){
    anchorColumnExists <- anchor %in% colnames(S4Vectors::mcols(gr))
    if (!anchorColumnExists){
        stop(sprintf("anchor '%s' is not found in mcols(gr)", anchor))
    } 
    return(anchor)
}


.validateTssWindow <- function(tss_window=NULL){
    if (is.null(tss_window)){
        tss_window <- .default_tss_window
    } else {
        hasValidLength <- length(tss_window) == 2
        isNumericVector <- is.vector(tss_window, mode="numeric")
        isAllIntegers <- all(tss_window == round(tss_window))
        isOrdered <- tss_window[[1]] <= tss_window[[2]]
        stopifnot("Invalid values for tss_window. See documentation." = {
            all(c(hasValidLength, isNumericVector, isAllIntegers, isOrdered))
        })
    }
    return(tss_window)
}


#' @importFrom BiocGenerics width colnames
#' @importFrom S4Vectors mcols
.validateTssObject <- function(tssObject
){
    if (!is.null(tssObject)){
        stopifnot("tssObject must be a GRanges object" = {
            .isGRanges(tssObject)
        })
        widths <- BiocGenerics::width(tssObject)
        stopifnot("All ranges in tssObject must have a width of 1" = {
            all(widths == 1)
        })
        mcolnames <- BiocGenerics::colnames(S4Vectors::mcols(tssObject))
        stopifnot("mcols(tssObject) must have an 'ID' column" = {
            "ID" %in% mcolnames
        })
    }   
    return(tssObject)
}


#' @importFrom S4Vectors isTRUEorFALSE
.checkBoolean <- function(argument,
                          value
){
    if (!S4Vectors::isTRUEorFALSE(value)){
        stop(sprintf("%s argument must be TRUE or FALSE", argument))
    }
    invisible(NULL)
}


.checkSingleInteger <- function(argument,
                                value,
                                null_ok=TRUE,
                                sign="any"
){
    isNull <- is.null(value) && null_ok
    isSingleIntValue <- is.vector(value, mode="numeric") &&
        length(value) == 1 &&
        value == round(value)
    hasValidSign <- switch(sign,
                           "positive" = value > 0,
                           "non-negative" = value >= 0,
                           "negative" = value < 0,
                           "non-positive" = value <= 0,
                           TRUE)
    if (!isNull && (!isSingleIntValue || !hasValidSign)){
        errorMessage <- paste("%s argument must be a single",
                              ifelse(sign == "any", "", sign),
                              "integer value or NULL")
        stop(sprintf(errorMessage, argument))
    }
    invisible(NULL)
}


.checkString <- function(argument,
                         value
){
    isCharacterVector <- is.vector(value, mode="character")
    hasLengthOne <- length(value) == 1
    if (!isCharacterVector || !hasLengthOne){
        stop(sprintf("%s argument must be a character string", argument))
    }
}


#' @importFrom methods is
#' @importFrom Biostrings DNA_BASES DNA_ALPHABET
.validateDNACharacterVariable <- function(seq,
                                          argument,
                                          len=NULL,
                                          nullOk=TRUE,
                                          exactBases=TRUE
){
    if (is.null(seq)){
        if (!nullOk){
            stop(sprintf("%s cannot be NULL", argument))
        } else {
            return("")
        }
    }
    if (methods::is(seq, "DNAString") && (is.null(len) || len == 1)){
        seq <- as.character(seq)
        return(seq)
    }
    if (is(seq, "DNAStringSet") && (is.null(len) || length(seq) == len)){
        seq <- as.character(seq)
    }
    seq <- toupper(seq)
    seq <- chartr("U", "T", seq)
    if (exactBases){
        pattern <- Biostrings::DNA_BASES
    } else {
        pattern <- grep('[A-Z]', Biostrings::DNA_ALPHABET, value=TRUE)
    }
    pattern <- paste0("^[", paste(pattern, collapse=""), "]+$")
    isNotCharacterVector <- !is.vector(seq, mode="character")
    hasBadSymbols <- !all(grepl(pattern, seq)) && all(nchar(seq) > 0)
    hasBadLength <- !is.null(len) && length(seq) != len
    if(anyNA(seq) || hasBadSymbols || isNotCharacterVector || hasBadLength){
        object <- ifelse(is.null(len), "vector", "string")
        stop(sprintf("%s is not a valid DNA character %s", argument, object))
    }
    return(seq)
}




# convert values in scientific notation to integers
.makeLongInteger <- function(n){
    as.integer(format(n, scientific=FALSE))
}



# Those functions are faster than
# their Biostrings counterpart
# for short sequences; useful for alignment.

.complement <- function(x){
    ## paste(grep("[A-Z]", Biostrings::DNA_ALPHABET, value=TRUE), collapse="")
    chartr("ACGTMRWSYKVHDBN", "TGCAKYWSRMBDHVN", x)
}
.reverse <- function(x){
    x <- strsplit(x, split="")
    x <- lapply(x, rev)
    x <- vapply(x, paste, collapse="", FUN.VALUE=character(1))
    return(x)
}
.revComp <- function(x){
    x <- .reverse(x) 
    x <- .complement(x)
    return(x)
}


#' @importFrom Biostrings DNAStringSet reverseComplement
.revCompBs <- function(x){
    x <- Biostrings::DNAStringSet(x)
    x <- Biostrings::reverseComplement(x)
    x <- as.character(x)
    return(x)
}


## not currently used
# adds/drops "chr" prefix in seqlevels, as required format differs by package
#' @importFrom GenomeInfoDb seqlevels renameSeqlevels
.toggleSeqlevels <- function(gr, dropChr=TRUE){
    stopifnot(.isGRanges(gr))
    stopifnot(is.logical(dropChr))
    if (dropChr){
        new_seqlevels <- gsub('chr', '', GenomeInfoDb::seqlevels(gr))
    } else {
        new_seqlevels <- paste0('chr', GenomeInfoDb::seqlevels(gr))
    }
    gr <- GenomeInfoDb::renameSeqlevels(gr, new_seqlevels)
    return(gr)
}


#' @importFrom crisprBase pamSide cutSites spacerLength pams weights motifs
.identicalNucleases <- function(nuc1,
                                nuc2,
                                checkSpacerLength=TRUE,
                                checkWeights=TRUE
){
    identicalPamSide <- identical(crisprBase::pamSide(nuc1),
                                  crisprBase::pamSide(nuc2))
    identicalCutSites <- identical(crisprBase::cutSites(nuc1),
                                   crisprBase::cutSites(nuc2))
    identicalSpacerLength <- identical(crisprBase::spacerLength(nuc1),
                                       crisprBase::spacerLength(nuc2))
    identicalPams <- identical(crisprBase::pams(nuc1, as.character=TRUE),
                               crisprBase::pams(nuc2, as.character=TRUE))
    identicalWeights <- identical(crisprBase::weights(nuc1),
                                  crisprBase::weights(nuc2))
    identicalMotifs <- identical(crisprBase::motifs(nuc1, as.character=TRUE),
                                 crisprBase::motifs(nuc2, as.character=TRUE))
    identicalNucleases <- identicalPamSide &&
        identicalCutSites &&
        identicalPams &&
        identicalMotifs 
    if (checkSpacerLength){
        identicalNucleases <- identicalNucleases && identicalSpacerLength
    }
    if (checkWeights){
        identicalNucleases <- identicalNucleases && identicalWeights
    }                
    return(identicalNucleases)
}  

# Remove NULLs from a list
compact <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}




