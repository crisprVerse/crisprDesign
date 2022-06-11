#' @title Restriction enzyme recognition sites in spacer sequences
#' 
#' @description Functions for identifying spacers sequences that contain
#'     restriction recognition sites for specific restriction enzymes.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param enzymeNames Character vector of enzyme names. 
#' @param patterns Optional named character vector for custom restriction site
#'     patterns. Vector names are treated as enzymes names. See example.
#' @param includeDefault Should commonly-used enzymes be included?
#'     TRUE by default.
#' @param flanking5,flanking3 Character string indicating the 5' or 3' flanking
#'     sequence, respectively, of the spacer sequence in the lentivial vector.
#' 
#' @return \code{\link{getRestrictionEnzymes}} returns a DataFrame indicating
#'     whether cutting sites for the specified enzymes are found in the gRNA
#'     cassette (flanking sequences + spacer sequences).
#'     
#'     \code{\link{addRestrictionEnzymes}} adds the annotation from
#'     \code{getRestrictionEnzymes} to \code{guideSet}.
#' 
#' @details Restriction enzymes are often used for cloning purpose during the
#'     oligonucleotide synthesis of gRNA lentiviral constructs. Consequently,
#'     it is often necessary to avoid restriction sites of the used restriction
#'     enzymes in and around the spacer sequences. The functions
#'     \code{addRestrictionEnzymes} and \code{getRestrictionEnzymes} allows for
#'     flagging problematic spacer sequences by searching for restriction sites
#'     in the [flanking5][spacer][flanking3] sequence.
#'     
#'     The following enzymes are included when \code{includeDefault=TRUE}:
#'     EcoRI, KpnI, BsmBI, BsaI, BbsI, PacI, and MluI.
#'     
#'     Custom recognition sequences in \code{patterns} may use the IUPAC
#'     nucleotide code, excluding symbols indicating gaps. Avoid providing
#'     enzyme names in \code{patterns} that are already included by default (if
#'     \code{includeDefault=TRUE}) or given by \code{enzymeNames}. Patterns
#'     with duplicated enzyme names will be silently ignored, even if the
#'     recognition sequence differs. See example.
#' 
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{\link{enzymeAnnotation}} to retrieve existing enzyme
#'     annotation from a \linkS4class{GuideSet} object.
#' 
#' @examples
#' data(SpCas9, package="crisprBase")
#' seq <- c("ATTTCCGGAGGCGAATTCGGCGGGAGGAGGAAGACCGG")
#' guideSet <- findSpacers(seq, crisprNuclease=SpCas9)
#' 
#' # Using default enzymes:
#' guideSet <- addRestrictionEnzymes(guideSet)
#' 
#' # Using custom enzymes:
#' guideSet <- addRestrictionEnzymes(guideSet,
#'                                   patterns=c(enz1="GGTCCAA",
#'                                              enz2="GGTCG"))
#' 
#' # Avoid duplicate enzyme names
#' guideSet <- addRestrictionEnzymes(guideSet,
#'                                   patterns=c(EcoRI="GANNTC")) # ignored
#' 
#' @export
addRestrictionEnzymes <- function(guideSet,
                                  enzymeNames=NULL,
                                  patterns=NULL,
                                  includeDefault=TRUE,
                                  flanking5="ACCG",
                                  flanking3="GTTT"
){
    guideSet <- .validateGuideSetOrPairedGuideSet(guideSet)
    if (.isGuideSet(guideSet)){
        out <- addRestrictionEnzymes_guideset(guideSet,
                                              enzymeNames=enzymeNames,
                                              patterns=patterns,
                                              includeDefault=includeDefault,
                                              flanking5=flanking5,
                                              flanking3=flanking3)
    } else if (.isPairedGuideSet(guideSet)){
        unifiedGuideSet <- .pairedGuideSet2GuideSet(guideSet)
        unifiedGuideSet <- addRestrictionEnzymes_guideset(unifiedGuideSet,
                                                          enzymeNames=enzymeNames,
                                                          patterns=patterns,
                                                          includeDefault=includeDefault,
                                                          flanking5=flanking5,
                                                          flanking3=flanking3)
        out <- .addColumnsFromUnifiedGuideSet(guideSet,
                                              unifiedGuideSet)
    }
    return(out)
}







#' @importFrom S4Vectors split mcols<-
addRestrictionEnzymes_guideset <- function(guideSet,
                                           enzymeNames=NULL,
                                           patterns=NULL,
                                           includeDefault=TRUE,
                                           flanking5="ACCG",
                                           flanking3="GTTT"
){
    enzymeAnnotation <- getRestrictionEnzymes(guideSet,
                                              enzymeNames=enzymeNames,
                                              patterns=patterns,
                                              includeDefault=includeDefault,
                                              flanking5=flanking5,
                                              flanking3=flanking3)
    dfs <- S4Vectors::split(enzymeAnnotation,
                            f=factor(rownames(enzymeAnnotation),
                                     levels=names(guideSet)))
    S4Vectors::mcols(guideSet)[["enzymeAnnotation"]] <- dfs
    return(guideSet)
}








#' @rdname addRestrictionEnzymes
#' @export
#' @importFrom S4Vectors DataFrame
getRestrictionEnzymes <- function(guideSet,
                                  enzymeNames=NULL,
                                  patterns=NULL,
                                  includeDefault=TRUE,
                                  flanking5="ACCG",
                                  flanking3="GTTT"
){
    guideSet <- .validateGuideSet(guideSet)
    enzymeMotifs <- .enzymeMotifs(includeDefault=includeDefault,
                                  enzymeNames=enzymeNames,
                                  patterns=patterns)
    spacers <- .spacersWithFlankingRegions(guideSet=guideSet,
                                           flanking5=flanking5,
                                           flanking3=flanking3)
    enzymeAnnotation <- lapply(enzymeMotifs,
                               grepl,
                               x=spacers)
    enzymeAnnotation <- DataFrame(enzymeAnnotation,
                                  row.names=names(guideSet))
    return(enzymeAnnotation)
}




.enzymeMotifs <- function(includeDefault,
                          enzymeNames,
                          patterns
){
    stopifnot("includeDefault must be TRUE or FALSE" = {
        includeDefault %in% c(TRUE, FALSE) && length(includeDefault) == 1
    })
    if (includeDefault){
        enzymeNames <- unique(c(.defaultEnzymeNames, enzymeNames))
    }
    motifs <- .getEnzymeMotifs(enzymeNames)
    motifs <- .addCustomEnzymes(motifs, patterns)
    stopifnot("no restriction enzymes found" = {
        length(motifs) > 0
    })
    motifs <- vapply(motifs, .enzymeMotif2RegexPattern, FUN.VALUE=character(1))
    return(motifs)
}


.defaultEnzymeNames <- c("EcoRI", "KpnI", "BsmBI",
                         "BsaI", "BbsI", "PacI","MluI")


#' @importFrom crisprBase motifs
.getEnzymeMotifs <- function(enzymeNames){
    data("restrictionEnzymes",
         package="crisprBase",
         envir=environment())
    .checkEnzymeNames(enzymeNames, restrictionEnzymes)
    motifs <- vapply(enzymeNames, function(x){
        crisprBase::motifs(restrictionEnzymes[[x]],
                           as.character=TRUE)
    }, FUN.VALUE=character(1))
    return(motifs)
}


.checkEnzymeNames <- function(enzymeNames,
                              restrictionEnzymes){
    if (length(enzymeNames) > 0){
        stopifnot("enzymeNames must be a character vector" = {
            is.vector(enzymeNames, mode="character")
        })
        badNames <- setdiff(enzymeNames, names(restrictionEnzymes))
        if (length(badNames) > 0){
            stop("restriction enzyme name(s) not found: ",
                 paste(badNames, collapse=', '))
        }
    }
    invisible(NULL)
}


.addCustomEnzymes <- function(motifs,
                              patterns
){
    if (!is.null(patterns)){
        patterns <- .validateCustomEnzymes(patterns)
        patterns <- as.list(patterns)
        patterns <- patterns[setdiff(names(patterns), names(motifs))]
        motifs <- c(motifs, patterns)
    }
    return(motifs)
}


.validateCustomEnzymes <- function(patterns
){
    stopifnot("patterns must be a character vector" = {
        is.vector(patterns, mode="character")
    })
    stopifnot("patterns vector must have names" = {
        !is.null(names(patterns)) &&
            !any(c(NA, "") %in% names(patterns))
    })
    stopifnot("patterns vector must have unique names" = {
        all.equal(names(patterns), unique(names(patterns)))
    })
    patterns <- .validateDNACharacterVariable(patterns,
                                              "patterns",
                                              exactBases=FALSE)
    return(patterns)
}



.enzymeMotif2RegexPattern <- function(motif
){
    revMotif <- .revCompBs(motif)
    pattern <- c(.iupacCode2RegexPattern(motif),
                 .iupacCode2RegexPattern(revMotif))
    pattern <- unique(pattern)
    pattern <- paste0(pattern, collapse="|")
    return(pattern)
}


#' @importFrom Biostrings DNA_BASES IUPAC_CODE_MAP
.iupacCode2RegexPattern <- function(seq
){
    seqBases <- strsplit(seq, '')[[1]]
    patternBases <- vapply(seqBases, function(x){
        if (!x %in% Biostrings::DNA_BASES){
            x <- paste0("[", Biostrings::IUPAC_CODE_MAP[x], "]")
        }
        x
    }, FUN.VALUE=character(1))
    pattern <- paste0(patternBases, collapse="")
    return(pattern)
}




.spacersWithFlankingRegions <- function(guideSet,
                                        flanking5,
                                        flanking3
){
    spacers    <- spacers(guideSet,
                          as.character=TRUE)
    flanking5 <- .validateDNACharacterVariable(flanking5, "flanking5", len=1)
    flanking3 <- .validateDNACharacterVariable(flanking3, "flanking3", len=1)
    spacers    <- paste0(flanking5, spacers, flanking3, recycle0=TRUE)
    names(spacers) <- names(guideSet)
    return(spacers)
}
