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
})
