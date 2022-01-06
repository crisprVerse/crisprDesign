# VCF_PATH <- "~/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"
VCF_PATH <- system.file("extdata",
                        file="common_snps_dbsnp151_example.vcf.gz",
                        package="crisprDesignS4")


test_that("guideSet argument must be a GuideSet object", {
    gs_as_gr <- GRanges(seqnames=seqnames(guideSetExample),
                        ranges=IRanges(start=start(guideSetExample),
                                       width=width(guideSetExample)),
                        strand=strand(guideSetExample))
    mcols(gs_as_gr) <- mcols(guideSetExample)
    names(gs_as_gr) <- names(guideSetExample)
    metadata(gs_as_gr) <- metadata(guideSetExample)
    bad_input <- list("guideSetExample",
                      as.data.frame(guideSetExample),
                      gs_as_gr)
    lapply(bad_input, function(x){
        expect_error(addSNPAnnotation(x, vcf=VCF_PATH))
    })
    expect_error(addSNPAnnotation(guideSetExample[1], vcf=VCF_PATH),
                 regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    expect_error(addSNPAnnotation(guideSetExample[0], vcf=VCF_PATH),
                 regexp=NA)
})


test_that("maf argument must be a numeric value between 0-1 (exclusive)", {
    skip("long run time")
    
    bad_values <- list(NULL,
                       NA,
                       "0.1",
                       -0.1,
                       1,
                       1.1,
                       c(0.1, 0.2))
    lapply(bad_values, function(x){
        expect_error(addSNPAnnotation(guideSetExample, vcf=VCF_PATH, maf=x))
    })
    good_values <- list(0,
                        0.1)
    lapply(good_values, function(x){
        expect_error(addSNPAnnotation(guideSetExample, vcf=VCF_PATH, maf=x),
                     regexp=NA)
    })
})


test_that("vcf argument must be a VCF object or a path to a vcf file", {
    bad_values <- list(NULL,
                       NA,
                       0,
                       "?badFilePath",
                       "bad/file.path",
                       c(VCF_PATH, VCF_PATH))
    lapply(bad_values, function(x){
        expect_error(addSNPAnnotation(guideSetExample, vcf=x))
    })
    expect_error(addSNPAnnotation(guideSetExample[1], vcf=VCF_PATH),
                 regexp=NA)
})


test_that("hasSNP logical column is appended to guideSet", {
    skip("long run time")
    
    out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH)
    expect_true("hasSNP" %in% colnames(mcols(out)))
    expect_true(is.logical(mcols(out)[["hasSNP"]]))
})


test_that("hasSNP value matches snps annotation list-column", {
    skip("long run time")
    
    out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH)
    hasSNP <- out$hasSNP
    
    snpAnnotation <- out[hasSNP]
    snpAnnotation <- vapply(snpAnnotation$snps, nrow, FUN.VALUE=numeric(1))
    expect_true(all(snpAnnotation > 0))
    
    noSnpAnnotation <- out[!hasSNP]
    noSnpAnnotation <- vapply(noSnpAnnotation$snps, nrow, FUN.VALUE=numeric(1))
    expect_true(all(noSnpAnnotation == 0))
})


test_that("SNP annotation (snps) as SplitDataFrameList is added to guideSet", {
    out <- addSNPAnnotation(guideSetExample[1], vcf=VCF_PATH)
    expect_true("snps" %in% colnames(mcols(out)))
    expect_true(is(mcols(out)[["snps"]], "SplitDataFrameList"))
})


test_that("guideSet with no SNP annotation is handled gracefully", {
    out <- addSNPAnnotation(guideSetExample[1], vcf=VCF_PATH)
    expect_true(nrow(snps(out)) == 0)
    expect_false(any(out$hasSNP))
})


test_that("SNP annotation is included for each spacer-SNP overlap", {
    guideSet <- guideSetExample[c(10, 10)]
    out <- addSNPAnnotation(guideSet, vcf=VCF_PATH)
    expect_identical(out$snps[[1]], out$snps[[2]])
    
    names(guideSet)[1] <- paste0(names(guideSet)[1], "_2")
    guideSet[1] <- GenomicRanges::shift(guideSet[1], 1)
    guideSet$pam_site[1] <- guideSet$pam_site[1] + 1
    guideSet$cut_site[1] <- guideSet$cut_site[1] + 1
    out <- addSNPAnnotation(guideSet, vcf=VCF_PATH)
    expect_identical(out$snps[[1]]$rs, out$snps[[2]]$rs)
    expect_false(identical(out$snps[[1]]$rs_site_rel, out$snps[[2]]$rs_site_rel))
})


test_that("SNP annotation only includes SNPs with MAF >= maf argument", {
    skip("long run time")
    
    maf <- list(0.01,
                0.5,
                0.9999)
    lapply(maf, function(x){
        out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH, maf=x)
        expect_true(nrow(snps(out)) == 0 ||
                        all(snps(out)$MAF_1000G >= x) ||
                        all(snps(out)$MAF_TOPMED >= x))
    })
})


test_that("SNP annotation correctly identifies 'snp' variants", {
    out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH)
    out <- snps(out)
    out <- out[out$type == "snp", , drop=FALSE]
    expect_true(all(out$length) == 1)
    expect_true(all(nchar(out$allele_ref) == nchar(out$allele_minor)))
})


test_that("SNP annotation correctly identifies 'del' variants", {
    skip("no del for guideSetExample")
    
    out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH)
    out <- snps(out)
    out <- out[out$type == "del", , drop=FALSE]
    expect_true(all(out$length) >= 1)
    expect_true(all(nchar(out$allele_ref) > nchar(out$allele_minor)))
})


test_that("SNP annotation correctly identifies 'ins' variants", {
    out <- addSNPAnnotation(guideSetExample, vcf=VCF_PATH)
    out <- snps(out)
    out <- out[out$type == "ins", , drop=FALSE]
    expect_true(all(out$length) >= 1)
    expect_true(all(nchar(out$allele_ref) < nchar(out$allele_minor)))
})
