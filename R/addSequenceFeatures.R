#' @title Add spacer sequence feature annotation columns to a
#'     \linkS4class{GuideSet} object
#' @description Add spacer sequence feature annotation columns, such as
#'     GC content, homopolymers, and hairpin predictions, to a
#'     \linkS4class{GuideSet} object.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param addHairpin Should predicted hairpin formation via
#'     sequence complementarity be calculated? FALSE by default.
#' @param backbone Backbone sequence in the guide RNA that is susceptible
#'     to hairpin formation with a complementary region in the spacer sequence
#'     (>= 5bp).
#' 
#' @return \code{guideSet} with the following columns appended to
#'    \code{mcols(guideSet)}: 
#'    \itemize{
#'        \item \code{percentGC} — percent GC content
#'        \item \code{polyA}, \code{polyC}, \code{polyG}, \code{polyT} —
#'        presence of homopolymers of 4nt or longer
#'        \item \code{selfHairpin} — prediction of hairpin formation within the
#'        spacer sequence via self-complementarity (>= 4bp)
#'        \item \code{backboneHairpin} — prediction of hairpin formation with
#'        the backbone sequence via complementarity (>= 5bp)
#'    }
#'
#' @examples
#' custom_seq <- c("ATTTCCGGAGGCGGAGAGGCGGGAGGAGCG")
#' data(SpCas9, package="crisprBase")
#' guideSet <- findSpacers(custom_seq, crisprNuclease=SpCas9)
#' guideSet <- addSequenceFeatures(guideSet)
#' 
#' 
#' 
#' @export
#' @importFrom S4Vectors mcols<-
addSequenceFeatures <- function(guideSet,
                                addHairpin=FALSE,
                                backbone="AGGCTAGTCCGT"
){
    guideSet <- .validateGuideSet(guideSet)
    
    seqs <- spacers(guideSet)
    guideSet <- .addPercentGC(guideSet, seqs)
    guideSet <- .addHomopolymers(guideSet, seqs)
    S4Vectors::mcols(guideSet)[["startingGGGGG"]] <- grepl("^GGGGG", seqs)
    
    stopifnot("addHairpin must be TRUE or FALSE" = {
        isTRUEorFALSE(addHairpin)
    })
    if (addHairpin){
        guideSet <- .addHairpins(guideSet, seqs, backbone)
    }
    
    return(guideSet)
}


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
.addHairpins <- function(guideSet,
                         seqs,
                         backbone
){
    # preG is optional (Thyme, 2016)
    seqs <- Biostrings::DNAStringSet(paste0("G", seqs, "GT")) 
    gap <- paste0(rep("-", spacerLength(guideSet) + 3), collapse="")
    backbone <- .validateDNACharacterVariable(seq=backbone,
                                              argument="backbone",
                                              len=1,
                                              nullOk=FALSE)
    backbone <- Biostrings::DNAString(paste0(gap, backbone))
    
    hairpins <- vapply(seqs, function(x){
        c(selfHairpin = .selfHairpin(x),
          backboneHairpin = .backboneHairpin(x, backbone))
    }, FUN.VALUE = logical(2))
    
    guideSet$selfHairpin <- hairpins['selfHairpin',]
    guideSet$backboneHairpin <- hairpins['backboneHairpin',]
    return(guideSet)
}


.selfHairpin <- function(seq
){
    min.armlength <- 4
    min.looplength <- 4
    max.looplength <- nchar(seq)
    hasSelfHairpin <- .findComplementarity(seq=seq,
                                           min.armlength=min.armlength,
                                           min.looplength=min.looplength,
                                           max.looplength=max.looplength)
    return(hasSelfHairpin)
}


#' @importFrom Biostrings DNAString
.backboneHairpin <- function(seq,
                             backbone
){
    min.armlength <- 5
    min.looplength <- nchar(seq)
    max.looplength <- nchar(seq) + nchar(backbone)
    seq <- Biostrings::DNAString(paste0(seq, backbone))
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
