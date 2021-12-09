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
    backbone <- .validateDNACharacterVariable(seq=backbone,
                                              argument="backbone",
                                              len=1,
                                              nullOk=FALSE)
    backbone <- Biostrings::DNAStringSet(backbone)
    
    hairpins <- vapply(seqs, function(x){
        c(selfHairpin = .selfHairpin(x),
          backboneHairpin = .backboneHairpin(x, backbone))
    }, FUN.VALUE = logical(2))
    
    guideSet$selfHairpin <- hairpins['selfHairpin',]
    guideSet$backboneHairpin <- hairpins['backboneHairpin',]
    return(guideSet)
}


#' @importFrom Biostrings nchar subseq reverseComplement
.selfHairpin <- function(seq
){
    minStemLength <- 4
    minLoopLength <- 4
    spacerWindowCount <- Biostrings::nchar(seq)-2*minStemLength-minLoopLength+1
    startingHairpinIndices <- seq_len(spacerWindowCount)
    for (i in startingHairpinIndices){
        queryStem <- Biostrings::subseq(seq, i, width=minStemLength)
        stemPercentGC <- .calculatePercentGC(queryStem)
        if (stemPercentGC >= 50){
            subjectStem <- Biostrings::subseq(seq,
                                              i+minStemLength+minLoopLength,
                                              Biostrings::nchar(seq))
            subjectStem <- Biostrings::reverseComplement(subjectStem)
            if (grepl(queryStem, subjectStem)){
                return(TRUE)
            }
        }
    }
    return(FALSE)
}


#' @importFrom Biostrings nchar reverseComplement subseq
.backboneHairpin <- function(seq,
                             backbone
){
    backbone <- Biostrings::reverseComplement(backbone)
    minStemLength <- 5
    spacerWindowCount <- Biostrings::nchar(seq) - minStemLength + 1
    startingHairpinIndices <- seq_len(spacerWindowCount)
    for (i in startingHairpinIndices){
        queryStem <- Biostrings::subseq(seq, i, width=minStemLength)
        stemPercentGC <- .calculatePercentGC(queryStem)
        if (stemPercentGC >= 50 && grepl(queryStem, backbone)){
            return(TRUE)
        }
    }
    return(FALSE)
}
