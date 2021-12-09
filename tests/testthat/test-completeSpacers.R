#' @param start Coordinate of the first nucleotide of the spacer sequences.
#'     Must be always less than \code{end}.
#' @param end Coordinate of the last nucleotide of the spacer sequence.
#'     Must be always greater than \code{start}.
#' @param chr The chromosome in which the protospacer sequence is located.
#' @param pam_site Coordinate of the first nucleotide of the PAM sequence.
#' @param strand Either "+" or "-".
#' @param crisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param genome Either hg38 or mm10.
#' @param spacerLen Spacer sequence length.
#'     If NULL, the information is obtained from \code{crisprNuclease}.
#' @param cut_offset Distance in nucleotides between \code{pam_site}
#'     and \code{cut_site}. If NULL, the information is obtained from 
#'     code{crisprNuclease}.

# getPAMSiteFromStartAndEnd <- function(start=NULL,
#                                       end=NULL,
#                                       strand,
#                                       crisprNuclease=NULL,
#                                       spacerLen=NULL



## tests for getPAMSiteFromStartAndEnd()


test_that("start and end args must be positive integers or NULL", {
    
})


test_that("at least one of start/end arg must be provided", {
    
})


test_that("strand arg must be '+' or '-'", {
    
})


test_that("start, end, and strand args must have same length", {
    
})


test_that("crisprNuclease arg must be a CrisprNuclease object", {
    
})


test_that("spacerLen arg must be a single positive integer or NULL", {
    
})


# test requiring start/end/spacerLen be in agreement?


test_that("getPAMSiteFromStartAndEnd returns correct PAM site(s)", { # break into multiple tests varying 1 thing
    # tests with single start/end/strand values
    # tests with vector of start/end/strand values
    # tests with SpCas9
    # tests with AsCas12a
    # tests with +/- strands
})




# getCutSiteFromPamSite <- function(pam_site,
#                                   strand,
#                                   crisprNuclease=NULL,
#                                   cut_offset=NULL


## tests for getCutSiteFromPamSite()


test_that("pam_site arg must be vector of positive integer(s", {
    
})


test_that("strand arg must be '+' or '-'", {

})


test_that("crisprNuclease arg must be a CrisprNuclease object", { # add function name
    
})


test_that("cut_offset arg must be a single integer value or NULL", {
    
})


test_that("pam_site and strand args must have same lengths", {
    
})


test_that("getCutSiteFromPamSite returns correct cut site(s)", { # break into multiple tests varying 1 thing
    # tests with single pam_site/strand values
    # tests with vector of pam_site/strand values
    # tests with SpCas9
    # tests with AsCas12a
    # tests with +/- strands
})






# getPAMSequence <- function(chr,
#                            pam_site,
#                            strand, 
#                            crisprNuclease=NULL,
#                            genome=NULL


## tests for getCutSiteFromPamSite()


test_that("chr arg must be a character vector", {
    
})


test_that("pam_site arg must be a vector of positive integer(s)", {
    
})


test_that("strand arg must be '+' or '-'", {
    
})


test_that("crisprNuclease arg must be a CrisprNuclease object", { # add function name
    
})


test_that("genome arg must be a BSgenome object or a permitted string", {
    shorthand_genomes <- c('hg38', 'mm10')
    
})


test_that("chr arg must be in genome", {
    
})


test_that("chr, pam_site, strand args must have the same length", {
    # when unlisted...?
})


test_that("getCutSiteFromPamSite requires pam_site coordinate to exist", {
    
})


test_that("getCutSiteFromPamSite returns correct PAM sequences", { # break into multiple tests varying 1 thing
    # tests with single chr/pam_site/strand values
    # tests with vector of chr/pam_site/strand values
    # tests with SpCas9
    # tests with AsCas12a
    # tests with +/- strands
    # test with human/mouse/other genomes
})




# getSpacerSequence <- function(chr,
#                               pam_site,
#                               strand,
#                               crisprNuclease=NULL,
#                               genome=NULL,
#                               spacerLen=NULL


## tests for getSpacerSequence()


test_that("chr arg must be a character vector", {
    
})


test_that("pam_site arg must be a vector of positive integer(s)", {
    
})


test_that("strand arg must be '+' or '-'", {
    
})


test_that("crisprNuclease arg must be a CrisprNuclease object", { # add function name
    
})


test_that("genome arg must be a BSgenome object or a permitted string", {
    shorthand_genomes <- c('hg38', 'mm10')
    
})


test_that("spacerLen arg must be a single positive integer or NULL", {
    
})


test_that("chr arg must be in genome", {
    
})


test_that("chr, pam_site, strand args must have the same length", {
    # when unlisted...?
})


test_that("getSpacerSequence requires pam_site coordinate to exist", {
    
})


test_that("getSpacerSequence returns correct spacer sequences", { # break into multiple tests varying 1 thing
    # tests with single chr/pam_site/strand values
    # tests with vector of chr/pam_site/strand values
    # tests with SpCas9
    # tests with AsCas12a
    # tests with +/- strands
    # test with human/mouse/other genomes
    # test with different spacerLen
})