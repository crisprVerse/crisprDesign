test_that("guideSet argument is required to be a GuideSet object", {
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
        expect_error(addOnTargetScores(x))
    })
    ## long run time
    # expect_error(addOnTargetScores(guideSetExample), regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    # expect_error(addOnTargetScores(guideSetExample[0]), regexp=NA)
    # expect_warning(addOnTargetScores(guideSetExample[0])) # scoring method not recognized: crisprscan
})


test_that("enzyme argument is required to be a specific value", {
    expect_error(addOnTargetScores(guideSetExample, enzyme="BAD_VALUE"))
    ## long run time
    # expect_error(addOnTargetScores(guideSetExample, enzyme=NULL), regexp=NA)
    # expect_error(addOnTargetScores(guideSetExample, enzyme="WT"), regexp=NA)
    # expect_error(addOnTargetScores(guideSetExample, enzyme="ESP"), regexp=NA)
    # expect_error(addOnTargetScores(guideSetExample, enzyme="HF"), regexp=NA)
})


test_that("promoter argument is required to be a specific value", {
    expect_error(addOnTargetScores(guideSetExample, promoter="BAD_VALUE"))
    ## long run time
    # expect_error(addOnTargetScores(guideSetExample, promoter=NULL),
    #              regexp=NA)
    # expect_error(addOnTargetScores(guideSetExample, promoter="U6"),
    #              regexp=NA)
    # expect_error(addOnTargetScores(guideSetExample, promoter="T7"),
    #              regexp=NA)
})


test_that("methods argument is required to match specific values", {
    bad_input <- list(NULL,
                      "BAD_VALUE",
                      c("azimuth", "BAD_VALUE"))
    lapply(bad_input, function(x){
        expect_error(addOnTargetScores(guideSetExample[1], methods=x))
    })
    ## long run time
    # good_input <- list("azimuth",
    #                    c("azimuth", "azimuth"),
    #                    c("azimuth", "deephf", "ruleset1", "lindel",
    #                      "deepcpf1", "enpamgb"))
    # lapply(good_input, function(x){
    #     expect_error(addOnTargetScores(guideSetExample[1], methods=x),
    #                  regexp=NA)
    # })
})


test_that("enzyme argument only affects deephf scores", {
    ## long run time
    # methods <- c("azimuth", "ruleset1", "lindel", "deepcpf1", "enpamgb")
    # out_wt <- addOnTargetScores(guideSetExample, enzyme="WT", methods=methods)
    # out_esp <- addOnTargetScores(guideSetExample, enzyme="ESP", methods=methods)
    # out_hf <- addOnTargetScores(guideSetExample, enzyme="HF", methods=methods)
    # expect_identical(out_wt, out_esp)
    # expect_identical(out_wt, out_hf)
    # out_deephf_wt <- addOnTargetScores(guideSetExample, enzyme="WT",
    #                                    methods="deephf")$score_deephf
    # out_deephf_esp <- addOnTargetScores(guideSetExample, enzyme="ESP",
    #                                    methods="deephf")$score_deephf
    # out_deephf_hf <- addOnTargetScores(guideSetExample, enzyme="HF",
    #                                    methods="deephf")$score_deephf
    # expect_false(identical(out_deephf_wt, out_deephf_esp))
    # expect_false(identical(out_deephf_wt, out_deephf_hf))
    # expect_false(identical(out_deephf_esp, out_deephf_hf))
})


test_that("promoter argument only affects deephf scores", {
    ## long run time
    # methods <- c("azimuth", "ruleset1", "lindel", "deepcpf1", "enpamgb")
    # out_u6 <- addOnTargetScores(guideSetExample, promoter="U6", methods=methods)
    # out_t7 <- addOnTargetScores(guideSetExample, enzyme="HF", methods=methods)
    # expect_identical(out_u6, out_t7)
    # out_deephf_u6 <- addOnTargetScores(guideSetExample, promoter="U6",
    #                                    methods="deephf")$score_deephf
    # out_deephf_t7 <- addOnTargetScores(guideSetExample, promoter="T7",
    #                                     methods="deephf")$score_deephf
    # expect_false(identical(out_deephf_u6, out_deephf_t7))
})


test_that("appropriate methods are selected via guideSet crisprNuclease", {
    ## long run time
    # score_cols <- sort(c("score_azimuth", "score_deephf",
    #                      "score_ruleset1", "score_lindel"))
    # outSpcas9 <- addOnTargetScores(guideSetExample[1])
    # outSpcas9 <- sort(grep("^score_", colnames(mcols(outSpcas9)), value=TRUE))
    # expect_identical(score_cols, outSpcas9)
    # 
    # iqsec3 <- queryTxObject(grListExample,
    #                         "transcripts",
    #                         "gene_symbol",
    #                         "IQSEC3")
    # genome(iqsec3) <- "hg38"
    # guides_ascas12a <- findSpacers(iqsec3,
    #                                bsgenome=BSgenome.Hsapiens.UCSC.hg38,
    #                                crisprNuclease=AsCas12a)[1]
    # guides_enascas12a <- findSpacers(iqsec3,
    #                                  bsgenome=BSgenome.Hsapiens.UCSC.hg38,
    #                                  crisprNuclease=enAsCas12a)[1]
    # 
    # score_cols <- sort(c("score_deepcpf1"))
    # outAscas12a <- addOnTargetScores(guides_ascas12a)
    # outAscas12a <- sort(grep("^score_", colnames(mcols(outAscas12a)),
    #                          value=TRUE))
    # expect_identical(score_cols, outAscas12a)
    # 
    # score_cols <- sort(c("score_deepcpf1", "score_enpamgb"))
    # outEnAscas12a <- addOnTargetScores(guides_enascas12a)
    # outEnAscas12a <- sort(grep("^score_", colnames(mcols(outEnAscas12a)),
    #                            value=TRUE))
    # expect_identical(score_cols, outEnAscas12a)
    # 
    # expect_error(addOnTargetScores(guideSetExample,
    #                                methods=c("deepcpf1", "enpamgb")))
    # expect_error(addOnTargetScores(guides_ascas12a,
    #                                methods=c("azimuth", "deephf",
    #                                          "ruleset1", "lindel",
    #                                          "enpamgb")))
    # expect_error(addOnTargetScores(guides_enascas12a,
    #                                methods=c("azimuth", "deephf",
    #                                          "ruleset1", "lindel")))
})


test_that("scores for spacers with insufficient context are NA", {
    ## long run time
    # custom_seq <- "CAGTCAGTCAGTCAGTCAGTAGG"
    # guide <- findSpacers(custom_seq, crisprNuclease=SpCas9)
    # out_short <- addOnTargetScores(guide)
    # expect_equal(out_short$score_azimuth, NA)
    # expect_false(is.na(out_short$score_deephf))
    # expect_equal(out_short$score_lindel, NA)
    # expect_equal(out_short$score_ruleset1, NA)
})


test_that("scores requiring canonical PAM return NA for noncanonical PAM", {
    # scoring method not recognized: crisprscan
    # guides <- guideSetExample[1:2]
    # guides$pam <- DNAStringSet(c("AAG", "AAA"))
    # expect_warning(out_pam <- addOnTargetScores(guides))
    # expect_true(all(is.na(out_pam$score_azimuth)))
    # expect_true(all(is.na(out_pam$score_deephf)))
    # expect_true(all(is.na(out_pam$score_lindel)))
    # expect_true(all(is.na(out_pam$score_ruleset1)))
    # 
    # iqsec3 <- queryTxObject(grListExample,
    #                         "transcripts",
    #                         "gene_symbol",
    #                         "IQSEC3")
    # genome(iqsec3) <- "hg38"
    
    # guides_ascas12a <- findSpacers(iqsec3,
    #                                bsgenome=BSgenome.Hsapiens.UCSC.hg38,
    #                                crisprNuclease=AsCas12a)[1:2]
    # guides_ascas12a$pam <- DNAStringSet(c("GGGG", "GGGG"))
    # out <- addOnTargetScores(guides_ascas12a)
    # expect_true(all(is.na(out$score_deepcpf1)))   # still computes scores
    
    # guides_enascas12a <- findSpacers(iqsec3,
    #                                  bsgenome=BSgenome.Hsapiens.UCSC.hg38,
    #                                  crisprNuclease=enAsCas12a)[1:2]
    # guides_enascas12a$pam <- DNAStringSet(c("AACC", "GGGG"))
    # out <- addOnTargetScores(guides_enascas12a)
    # expect_true(all(!is.na(out$score_deepcpf1)))
    # expect_true(all(!is.na(out$score_enpamgb)))
})


test_that("spacers are annotated with correct score values (indexing)", {
    ## long run time
    # seqs <- paste0(spacers(guideSetExample), pams(guideSetExample))
    # scores <- getDeepHFScores(seqs)
    # out <- addOnTargetScores(guideSetExample, method="deephf")
    # expect_identical(scores$score, mcols(out)$score_deephf)
})

