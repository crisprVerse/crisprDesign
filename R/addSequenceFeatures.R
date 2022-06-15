#' @title Add spacer sequence feature annotation columns to a
#'     \linkS4class{GuideSet} object
#' @description Add spacer sequence feature annotation columns, such as
#'     GC content, homopolymers, and hairpin predictions, to a
#'     \linkS4class{GuideSet} object.
#' 
#' @param object A \linkS4class{GuideSet} or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param addHairpin Whether to include predicted hairpin formation via
#'     sequence complementarity. FALSE by default. See details.
#' @param backbone Backbone sequence in the guide RNA that is susceptible
#'     to hairpin formation with a complementary region in the spacer sequence.
#'     
#' @details \code{addSequenceFeatures} predicts spacers to form internal
#'     hairpins when there is a palindromic sequence within the spacer having
#'     arms of >=4nt and >=50\% GC content, and are separated by a loop of
#'     >=4nt. Backbone hairpin formation is predicted when the spacer and
#'     backbone share a complementary sequence of >=5nt and >=50\% GC content. 
#' 
#' @return \code{guideSet} with the following columns appended to
#'    \code{mcols(guideSet)}: 
#'    \itemize{
#'        \item \code{percentGC} — percent GC content
#'        \item \code{polyA}, \code{polyC}, \code{polyG}, \code{polyT} —
#'        presence of homopolymers of 4nt or longer
#'        \item \code{selfHairpin} — prediction of hairpin formation within the
#'        spacer sequence via self-complementarity
#'        \item \code{backboneHairpin} — prediction of hairpin formation with
#'        the backbone sequence via complementarity
#'    }
#'
#' @examples
#' custom_seq <- c("ATTTCCGGAGGCGGAGAGGCGGGAGGAGCG")
#' data(SpCas9, package="crisprBase")
#' guideSet <- findSpacers(custom_seq, crisprNuclease=SpCas9)
#' guideSet <- addSequenceFeatures(guideSet)
#' 
#' 
#' @rdname addSequenceFeatures
#' @importFrom S4Vectors mcols<-
setMethod("addSequenceFeatures", "GuideSet", function(object,
                                                      addHairpin=FALSE,
                                                      backbone="AGGCTAGTCCGT"
){
    object <- .validateGuideSet(object)
    
    seqs <- spacers(object)
    object <- .addPercentGC(object, seqs)
    object <- .addHomopolymers(object, seqs)
    S4Vectors::mcols(object)[["startingGGGGG"]] <- grepl("^GGGGG", seqs)
    
    stopifnot("addHairpin must be TRUE or FALSE" = {
        isTRUEorFALSE(addHairpin)
    })
    if (addHairpin){
        object <- .addHairpins(object, seqs, backbone)
    }
    
    return(object)
})


#' @rdname addSequenceFeatures
#' @export
setMethod("addSequenceFeatures",
          "PairedGuideSet", function(object,
                                     addHairpin=FALSE,
                                     backbone="AGGCTAGTCCGT"
){
    object <- .validateGuideSetOrPairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addSequenceFeatures(unifiedGuideSet,
                                           addHairpin=addHairpin,
                                           backbone=backbone)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})


#' @rdname addSequenceFeatures
#' @export
setMethod("addSequenceFeatures", "NULL", function(object){
    return(NULL)
})





#' @importFrom S4Vectors mcols<-
.addPercentGC <- function(guideSet,
                          seqs
){
    S4Vectors::mcols(guideSet)[['percentGC']] <- .calculatePercentGC(seqs)
    return(guideSet)
}


#' @importFrom Biostrings DNAStringSet letterFrequency
.calculatePercentGC <- function(seqs
){
    seqs <- Biostrings::DNAStringSet(seqs)
    percentGC <- Biostrings::letterFrequency(seqs, c("C", "G"), as.prob=TRUE)
    percentGC <- rowSums(percentGC)
    percentGC <- round(100*percentGC, 1)
    return(percentGC)
}


#' @importFrom Biostrings DNA_BASES
#' @importFrom S4Vectors mcols<-
.addHomopolymers <- function(guideSet,
                             seqs
){
    bases <- Biostrings::DNA_BASES
    for (i in bases){
        pattern <- paste0(rep(i, 4), collapse='')
        hasHomopolymer <- grepl(pattern, seqs, ignore.case=TRUE)
        colname <- paste0('poly', i)
        S4Vectors::mcols(guideSet)[[colname]] <- hasHomopolymer
    }
    return(guideSet)
}


#' @importFrom Biostrings DNAStringSet
#' @importFrom S4Vectors mcols<-
.addHairpins <- function(guideSet,
                         seqs,
                         backbone
){
    # preG is optional (Thyme, 2016)
    seqs <- paste0("G", seqs, "GT", recycle0=TRUE)
    seqs <- Biostrings::DNAStringSet(seqs) 
    backbone <- .validateDNACharacterVariable(seq=backbone,
                                              argument="backbone",
                                              len=1,
                                              nullOk=FALSE)
    
    selfHairpin <- vapply(seqs, .selfHairpin, FUN.VALUE=logical(1))
    S4Vectors::mcols(guideSet)[["selfHairpin"]] <- selfHairpin
    
    gapLength <- spacerLength(guideSet) + 3 + S4Vectors::nchar(backbone)
    gap <- strrep("-", gapLength)
    seqs <- paste0(seqs, gap, backbone, recycle0=TRUE)
    seqs <- Biostrings::DNAStringSet(seqs)
    backboneHairpin <- vapply(seqs, .backboneHairpin, FUN.VALUE=logical(1))
    S4Vectors::mcols(guideSet)[["backboneHairpin"]] <- backboneHairpin
    
    return(guideSet)
}


.selfHairpin <- function(seq
){
    min.armlength <- 4
    min.looplength <- 4
    max.looplength <- max(nchar(seq), min.looplength)
    hasSelfHairpin <- .findComplementarity(seq=seq,
                                           min.armlength=min.armlength,
                                           min.looplength=min.looplength,
                                           max.looplength=max.looplength)
    return(hasSelfHairpin)
}


#' @importFrom Biostrings letterFrequency
.backboneHairpin <- function(seq
){
    min.armlength <- 5
    min.looplength <- Biostrings::letterFrequency(seq, letters="-")
    max.looplength <- max(nchar(seq), min.looplength)
    hasBackboneHairpin <- .findComplementarity(seq=seq,
                                               min.armlength=min.armlength,
                                               min.looplength=min.looplength,
                                               max.looplength=max.looplength)
    return(hasBackboneHairpin)
}


#' @importFrom Biostrings findPalindromes palindromeLeftArm
.findComplementarity <- function(seq,
                                 min.armlength,
                                 min.looplength,
                                 max.looplength
){
    palindromes <- Biostrings::findPalindromes(seq,
                                               min.armlength=min.armlength,
                                               min.looplength=min.looplength,
                                               max.looplength=max.looplength,
                                               max.mismatch=0,
                                               allow.wobble=FALSE)
    stem <- Biostrings::palindromeLeftArm(palindromes)
    stem <- as.character(stem)
    hasComplementarity <- any(.calculatePercentGC(stem) >= 50)
    return(hasComplementarity)
}
