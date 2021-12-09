test_that("guideSet argument requires a GuideSet object", {
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
        expect_error(addPamScores(x))
    })
    expect_error(addPamScores(guideSetExample), regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    expect_error(addPamScores(guideSetExample[0]), regexp=NA)
})


test_that("PAM scores are appended to guideSet", {
    out <- addPamScores(guideSetExample[1])
    expect_true("score_pam" %in% colnames(mcols(out)))
})


test_that("correct PAM scores are returned from guideSet's crisprNuclease", {
    # test case where NA is returned (pam not documented in crisprNuclease object)
    # test with SpCas9, AsCas12a, enAsCas12a
})