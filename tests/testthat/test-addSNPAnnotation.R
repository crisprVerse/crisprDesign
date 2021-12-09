#' @param guideSet A \linkS4class{GuideSet} object.
#' @param maf Minimum minor allele frequency to report
#'     (for a least one source among 1000Genomes and TOPMED).
#'     Default value of 0.01.
#' @param vcf Either string specfying path to a VCF file,
#'     or a \linkS4class{VCF} object. If NULL, the function will try to find
#'     a default VCF file for common SNPS among files located at the
#'     \code{getCrisprIndicesDir} location. 
#' @param verbose Should messages be printed? 
#'     TRUE by default.











# addSNPAnnotation <- function(guideSet,
#                              maf=0.01,
#                              vcf=NULL,
#                              verbose=TRUE



test_that("guideSet argument must be a GuideSet object", {
    
})


test_that("maf argument must be a numeric value between 0-1", {
    
})


test_that("vcf argument must be a VCF object or a path to a vcf file", {
    
})


test_that("verbose argument must be a logical value", {
    
})


test_that("hasSNP logical column is appended to guideSet", {
    
})


test_that("hasSNP value matches snps annotation list-column", {
    
})


test_that("SNP annotation (snps) as SplitDataFrameList is added to guideSet", {
    
})


test_that("SNP annotation is included for each spacer-SNP overlap", {
    # i.e., SNP will appear twice in snps() if it overlaps two spacers in guideSet
})


test_that("SNP annotation only includes SNPs with MAF >= maf argument", {
    # applies to either of MAF_1000G or MAF_TOPMED 
})


test_that("SNP annotation correctly identifies 'snp' variants", {
    # check allele_ref, allele_minor, snp_type, snp_length
})


test_that("SNP annotation correctly identifies 'del' variants", {
    # check allele_ref, allele_minor, snp_type, snp_length
})


test_that("SNP annotation correctly identifies 'ins' variants", {
    # check allele_ref, allele_minor, snp_type, snp_length
})


test_that("verbose prints message when TRUE", {
    
})


