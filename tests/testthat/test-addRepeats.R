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
        expect_error(addRepeats(x, grRepeatsExample))
        expect_error(removeRepeats(x, grRepeatsExample))
    })
    expect_error(addRepeats(guideSetExample, grRepeatsExample), regexp=NA)
    expect_error(removeRepeats(guideSetExample, grRepeatsExample), regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    expect_error(addRepeats(guideSetExample[0], grRepeatsExample),
                 regexp=NA)
    expect_error(removeRepeats(guideSetExample[0], grRepeatsExample),
                 regexp=NA)
})


test_that("gr.repeats argument is required to be a GRanges object", {
    bad_input <- list("grRepeatsExample",
                      as.data.frame(grRepeatsExample),
                      list(grRepeatsExample))
    lapply(bad_input, function(x){
        expect_error(addRepeats(guideSetExample, x))
        expect_error(removeRepeats(guideSetExample, x))
    })
    expect_error(addRepeats(guideSetExample, grRepeatsExample), regexp=NA)
    expect_error(removeRepeats(guideSetExample, grRepeatsExample), regexp=NA)
})


test_that("an empty gr.repeats argument returns no overlaps with repeats", {
    inRepeats <- addRepeats(guideSetExample, grRepeatsExample[0])$inRepeats
    expect_false(any(inRepeats))
    expect_equal(length(guideSetExample),
                 length(removeRepeats(guideSetExample, grRepeatsExample[0])))
})


test_that("ignore.strand argument is required to be TRUE or FALSE", {
    good_input <- list(TRUE, FALSE)
    lapply(good_input, function(x){
        expect_error(addRepeats(guideSetExample, grRepeatsExample,
                                ignore.strand=x),
                     regexp=NA)
        expect_error(removeRepeats(guideSetExample, grRepeatsExample,
                                   ignore.strand=x),
                     regexp=NA)
    })
    bad_input <- list(NA, NULL, "TRUE", 0)
    lapply(bad_input, function(x){
        expect_error(addRepeats(guideSetExample, grRepeatsExample,
                                ignore.strand=x))
        expect_error(removeRepeats(guideSetExample, grRepeatsExample,
                                   ignore.strand=x))
    })
})


test_that("inRepeats column is appended to guideSet", {
    addRepeats_output <- addRepeats(guideSetExample, grRepeatsExample)
    removeRepeats_output <- removeRepeats(guideSetExample, grRepeatsExample)
    expect_true("inRepeats" %in% names(mcols(addRepeats_output)))
    expect_true("inRepeats" %in% names(mcols(removeRepeats_output)))
})


test_that("inRepeats correctly identifies spacers in repeat regions", {
    # compare output to known values or
    # compare output inRepeats to findOverlaps indices
})


test_that("ignore.strand is used as intended in determining repeat overlaps", {
    # guideSet and repeats on same strand (both) and opposite strands (both)
})


test_that("removeRepeats only retains spacers not in repeat regions", {
    # compare addRepeats and removeRepeats output (subset by names)
    # example that removes all spacers
})
