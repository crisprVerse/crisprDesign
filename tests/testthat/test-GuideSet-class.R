#' @param spacers Character vector of spacer sequences.
#' @param pams Character vector of PAM sequences.
#' @param seqnames Character vector of chromosome names.
#' @param pam_site Integer vector of PAM site coordinates.
#' @param strand Character vector of gRNA strand.
#'    Only accepted values are "+" and "-".
#' @param CrisprNuclease \linkS4class{CrisprNuclease} object.
#' @param genome String specifying genome (e.g, "mm10" or "hg38").
#' @param ... Additional arguments for class-specific methods


# GuideSet <- function(spacers = NA_character_,
#                      pams = NULL,
#                      seqnames = NA_character_,
#                      pam_site = 0L,
#                      strand = NA_character_,
#                      CrisprNuclease = NULL,
#                      genome=NA




# requires seqnames, spacers, CrisprNuclease, pams
# strand must be one of '+' '-' '*'
# pam_site can be negative...shouldn't allow? should a positive integer be required? (<2**31)
# genome is stored in metadata(GuideSet), not genome(GuideSet); no restrictions on value
# spacers, pams should be same length, with seqnames, strand, pam_site also same length, or length 1 (or NULL)



# gs <- GuideSet(seqnames='chr1', spacers=DNAStringSet('ACGTCATGACTCTACTACGT'), CrisprNuclease = SpCas9, pams=DNAStringSet(c('AGG', 'TGG')))





# methods:
# crisprNuclease(object)
# spacers(object, as.character = FALSE)
# pams(object, as.character = FALSE)
# pamSites(object)
# cutSites(object)
# protospacers(object, as.character = FALSE)
# spacerLength(object)
# protospacerLength(object)
# prototypeSequence(object)
# pamLength(object)
# pamSide(object)
# snps(object, unlist = TRUE)
# alignments(object, columnName="alignments", unlist=TRUE)
# onTargets(object, columnName="alignments", unlist=TRUE)
# offTargets(object, columnName="alignments", max_mismatches=Inf, unlist=TRUE)
# geneAnnotation(object, unlist=TRUE, gene_id=NULL, tx_id=NULL, gene_symbol=NULL)
# tssAnnotation(object, unlist=TRUE, gene_id=NULL, gene_symbol=NULL)
# enzymeAnnotation(object, unlist=TRUE)




test_that("object argument is required to have type 'GuideSet'", {
    gs_as_gr <- GRanges(seqnames=seqnames(guideSetExample),
                        ranges=IRanges(start=start(guideSetExample),
                                       width=width(guideSetExample)))
    mcols(gs_as_gr) <- mcols(guideSetExample)
    names(gs_as_gr) <- names(guideSetExample)
    
    accessors <- list(crisprNuclease,
                      spacers,
                      pams,
                      pamSites,
                      cutSites,
                      protospacers,
                      spacerLength,
                      protospacerLength,
                      prototypeSequence,
                      pamLength,
                      pamSide,
                      snps,
                      alignments,
                      onTargets,
                      offTargets,
                      geneAnnotation,
                      tssAnnotation,
                      enzymeAnnotation)
    
    lapply(accessors, function(fun){
        expect_error(fun("guideSetExample"))
        expect_error(fun(as.data.frame(guideSetExample)))
        expect_error(fun(gs_as_gr))
        expect_error(fun(guideSetExample), regexp=NA)
    })
})


test_that("annotation accessors return NULL when lacking that annotation", {
    accessors <- list(snps,
                      alignments,
                      onTargets,
                      offTargets,
                      geneAnnotation,
                      tssAnnotation,
                      enzymeAnnotation)
    lapply(accessors, function(fun){
        expect_null(fun(guideSetExample))
    })
})


# accessors that depend(?) on crisprNuclease ... compare with SpCas9 for example, make GuideSet using AsCas12a
# accessors: crisprNuclease, pamLength, pamSide, spacerLength



