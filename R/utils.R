#' @importFrom methods is
NULL

utils::globalVariables(c("restrictionEnzymes",
                         "scoringMethodsInfo"))
utils::globalVariables(c("SpCas9",
                         "AsCas12a",
                         "enAsCas12a"))



#' @export
#' @importFrom GenomeInfoDb seqnames
GenomeInfoDb::seqnames

#' @export
#' @importFrom S4Vectors mcols
S4Vectors::mcols


.default_tss_window <- c(-500, 500)

.isGRanges <- function(dat){
    return(is(dat, "GRanges"))
}

.isGRangesList <- function(x){
    is(x, "GRangesList")
}

.isTxDb <- function(x){
    is(x, "TxDb")
}




#' @importFrom methods is
.validateCrisprNuclease <- function(crisprNuclease){
    if (is.null(crisprNuclease)){
        crisprNuclease <- .getDefaultCrisprNuclease()
    } else {
        if (!is(crisprNuclease, "CrisprNuclease")){
            stop("Provided nuclease must be a 'CrisprNuclease' object. ")
        }
    }
    return(crisprNuclease)
}





#' @importFrom utils data
.getDefaultCrisprNuclease <- function(type=c("Cas9", "Cas12a")){
    type <- match.arg(type)
    if (type=="Cas9"){
        data("SpCas9",
             package="crisprBase",
             envir=environment())
        nuc <- SpCas9
    } else {
        data("AsCas12a",
             package="crisprBase",
             envir=environment())
        nuc <- AsCas12a
    }
    return(nuc)
}

#' @importFrom methods is
.validateBSgenome <- function(bsgenome){
    if (is.null(bsgenome)){
        bsgenome <- .getDefaultBSgenome()
    } else {
        if (!is(bsgenome, "BSgenome")){
            stop("Provided genome must be a 'BSgenome' object. ")
        }
    }
    return(bsgenome)
}

.getDefaultBSgenome <- function(){
    if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
        out <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 
    } else {
        out <- NULL
    }
    return(out)
}

.getBSGenome <- function(genome){
    choices <- c("hg38", "mm10")
    if (!genome %in% choices){
        stop("Could not automatically find a corresponding BSgenome object ",
             "for genome ", genome, ". Please specify BSgenome object using ",
             "the 'bsgenome' argument.")
    } else if (genome=='hg38'){
        if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38")){
            out <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        } else {
            stop("BSgenome.Hsapiens.UCSC.hg38 must be installed.")
        }
    } else if (genome=="mm10"){
        if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10")){
            out <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        } else {
            stop("BSgenome.Mmusculus.UCSC.mm10 must be installed.")
        }
    } 
    return(out)
}


.as_df <- function(dat){
  # check input
    if (!is.data.frame(dat) && !.isGRanges(dat)){
        return(dat)
    }
    df <- as.data.frame(dat,
                        stringsAsFactors=FALSE)
    df$chr <- df$seqnames
    df$seqnames <- NULL
    df$width <- NULL
    return(df)
}


#' @importFrom S4Vectors DataFrame
#' @importFrom BiocGenerics rownames
.asDataFrame <- function(data
){
    spacerIds <- names(data)
    # intermediate data.frame preserves GPos/GRanges formatting
    data <- as.data.frame(data)
    colnames(data)[colnames(data) == "seqnames"] <- "chr"
    colnames(data)[colnames(data) == "pos"] <- "anchor_site"
    data <- S4Vectors::DataFrame(data)
    BiocGenerics::rownames(data) <- spacerIds
    return(data)
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
.as_gr <- function(dat){
    # check input
    if (!is.data.frame(dat) && !.isGRanges(dat)){
        stop('Invalid input; dat must be either",
             "GRanges or data.frame.')
    }
    if (.isGRanges(dat)){
        return(dat)
    }
    if (!"start" %in% colnames(dat) &&
        !"pam_site" %in% colnames(dat)){
        stop("start or pam_site must be a column name.")
    } 
    if (!"start" %in% colnames(dat)){
        dat$start <- dat$end <- dat$pam_site
    }
    gr <- makeGRangesFromDataFrame(dat,
                                   keep.extra.columns=TRUE)
    return(gr)
}



.isGuideSet <- function(object){
    is(object, "GuideSet")
}


.validateGuideSet <- function(object){
    if (!.isGuideSet(object)){
        stop("Object must be a GuideSet")
    }
    return(object)
}


.validateMetCutoff <- function(met_cutoff){
    stopifnot(is.numeric(met_cutoff) && length(met_cutoff)==1)
    if (met_cutoff != Inf){
        if (met_cutoff < 0 || met_cutoff %% 1 != 0){
            stop('met_cutoff must be a positive integer.')
        }
    }
    return(met_cutoff)
}

.validateAnchor <- function(anchor, gr){
    if (!anchor %in% colnames(mcols(gr))){
        stop(sprintf("anchor '%s' is not found in mcols(gr).", anchor))
    } 
    return(anchor)
}



.validateTssWindow <- function(tss_window=NULL){
    if (is.null(tss_window)){
        tss_window <- .default_tss_window
    } else {
        validLength <- length(tss_window) == 2
        validType <- is.vector(tss_window, mode="numeric")
        validIntegers <- all(unlist(tss_window) %% 1 == 0)
        validOrder <- tss_window[[1]] <= tss_window[[2]]
        stopifnot("Invalid values for tss_window. See documentation." = {
            all(c(validLength, validType, validIntegers, validOrder))
        })
    }
    return(tss_window)
}


.validateGRanges <- function(obj){
    if (!is(obj, "GRanges")){
        stop("Object must be a GRanges. ")
    } 
    return(obj)
}


.validateGRangesList <- function(obj){
    if (.isTxDb(obj)){
        return(TxDb2GRangesList(obj))
    }
    stopifnot("Object must be a TxDb or GRangesList" = {
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



#' @importFrom BiocGenerics width colnames
#' @importFrom S4Vectors mcols
.validateTssObject <- function(tssObject
){
    if (!is.null(tssObject)){
        stopifnot("tssObject must be a GRanges object" = {
            is(tssObject, "GRanges")
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





.validateGRangesNames <- function(gr){
    grNames <- names(gr)
    suggestedNames <- paste0("region_",
                             seq_along(gr))
    if (is.null(grNames) | all(is.na(grNames)) | all(grNames=="")){
        names(gr) <- suggestedNames
    } else if (sum(duplicated(grNames))>0){
        stop("The GRanges object has duplicated names.")
    }
    return(gr)
}




.validateCustomSeqNames <- function(x){
    temp <- names(x)
    if (is.null(temp)){
        names(x) <- paste0("region_", seq_along(x))
    } else if (sum(duplicated(temp))>0){
        stop("The character vector has duplicated names.")
    }
    return(x)
}





.validateInputFindSpacers <- function(dna){
    if (!is(dna, "DNAStringSet")){
        stop("dna is not a DNAStringSet")
    }
    if (is.null(mcols(dna))){
        stop("mcols(dna) is NULL")
    }
    cols <- c("seqnames", "start", "end")
    if (!all(cols %in% colnames(mcols(dna)))){
        stop("seqnames, start, and end must be columns of mcols(dna)")
    }
    if (is.null(metadata(dna))){
        stop("metadata(dna) is NULL.")
    }
    if (!"genome" %in% names(metadata(dna))){
        stop("genome must be an element of metadata(dna)")
    }
    return(dna)
}



#' @importFrom S4Vectors isTRUEorFALSE
.checkSingleBoolean <- function(argument,
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
        value %% 1 == 0
    hasGoodSign <- switch(sign,
                          "positive" = value > 0,
                          "non-negative" = value >= 0,
                          "negative" = value < 0,
                          "non-positive" = value <= 0,
                          TRUE)
    if (!isNull && (!isSingleIntValue || !hasGoodSign)){
        errorMessage <- paste("%s argument must be a single",
                              ifelse(sign == "any", "", sign),
                              "integer value or NULL")
        stop(sprintf(errorMessage, argument))
    }
    invisible(NULL)
}




#' @importFrom Biostrings DNA_BASES DNA_ALPHABET
.validateDNACharacterVariable <- function(seq,
                                          argument,
                                          len=NULL,
                                          nullOk=TRUE,
                                          exactBases=TRUE
){
    if (nullOk && is.null(seq)){
        return("")
    }
    seq <- toupper(seq)
    seq <- chartr("U", "T", seq)
    if (exactBases){
        pattern <- Biostrings::DNA_BASES
    } else {
        pattern <- grep('[A-Z]', Biostrings::DNA_ALPHABET, value=TRUE)
    }
    pattern <- paste0("^[", paste(pattern, collapse=""), "]+$")
    hasBadSymbols <- !all(grepl(pattern, seq)) && all(nchar(seq) > 0)
    isNotCharacterVector <- !is.vector(seq, mode="character")
    hasBadLength <- !is.null(len) && length(seq) != len
    if(anyNA(seq) || hasBadSymbols || isNotCharacterVector || hasBadLength){
        object <- ifelse(is.null(len), "vector", "string")
        stop(argument, " is not a valid DNA character ", object)
    }
    return(seq)
}




#' @importFrom crisprBase pams
.isValidPAM <- function(chr,
                        pam_site,
                        strand,
                        genome=c("hg38", "mm10"),
                        crisprNuclease=NULL,
                        canonical=TRUE
){
    crisprNuclease <- .validateCrisprNuclease(crisprNuclease)
    pams <- getPAMSequence(chr=chr,
                           pam_site=pam_site,
                           strand=strand,
                           genome=genome,
                           crisprNuclease=crisprNuclease)
    choices <- pams(crisprNuclease, primary=canonical)
    wh <- pams %in% choices
    return(wh)
}




# Remove NULLs from a list
compact <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}






# convert values in scientific notation to integers
.makeLongInteger <- function(n){
    as.integer(format(n, scientific=FALSE))
}


.destring <- function(x){
    x <- strsplit(x, split="_")
    x <- lapply(x, function(y) y[[1]])
    unlist(x)
}

# functions to get reverse complement of DNA sequence
.complement <- function(x){
    chartr("ATGC","TACG",x)
}
.reverse <- function(x){
    x <- strsplit(x, split="")[[1]]
    x <- rev(x)
    x <- paste(x, collapse="")
    return(x)
}
.revComp <- function(x){
    x <- .reverse(x) 
    x <- .complement(x)
    return(x)
}

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



.identicalNucleases <- function(nuc1, nuc2){
    cond0 <- identical(pamSide(nuc1),
                       pamSide(nuc2))
    cond1 <- identical(cutSites(nuc1),
                       cutSites(nuc2))
    cond2 <- identical(spacerLength(nuc1),
                       spacerLength(nuc2))
    cond3 <- identical(pams(nuc1, as.character=TRUE),
                       pams(nuc2, as.character=TRUE))
    cond4 <- identical(weights(nuc1),
                       weights(nuc2))
    cond5 <- identical(motifs(nuc1, as.character=TRUE),
                       motifs(nuc2, as.character=TRUE))
    cond0 & cond1 & cond2 & cond3 & cond4 & cond5 
}  




