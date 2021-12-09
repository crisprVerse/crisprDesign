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
        expect_error(addTssAnnotation(x, tssObject=tssObjectExample))
    })
    expect_error(addTssAnnotation(guideSetExample, tssObject=tssObjectExample),
                 regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    expect_error(addTssAnnotation(guideSetExample[0],
                                  tssObject=tssObjectExample),
                 regexp=NA)
})


test_that("tssObject arg must be a GRanges object", {
    bad_input <- list(NULL,
                      grListExample,
                      list(tssObjectExample),
                      as.data.frame(tssObjectExample))
    lapply(bad_input, function(x){
        expect_error(addTssAnnotation(guideSetExample, tssObject=x))
    })
    expect_error(addTssAnnotation(guideSetExample, tssObject=tssObjectExample),
                 regexp=NA)
})


test_that("anchor arg must be either 'cut_site' or 'pam_site'", {
    expect_error(addTssAnnotation(guideSetExample,
                                  tssObject=tssObjectExample,
                                  anchor="BAD_VALUE"))
    expect_error(addTssAnnotation(guideSetExample,
                                  tssObject=tssObjectExample,
                                  anchor="cut_site"),
                 regexp=NA)
    expect_error(addTssAnnotation(guideSetExample,
                                  tssObject=tssObjectExample,
                                  anchor="pam_site"),
                 regexp=NA)
})


test_that("tss_window arg must be a length 2 numeric vector", {
    bad_input <- list(c("-500", "500"),
                      c(-500.5, 500.5),
                      c(-500),
                      c(-500, 0, 500),
                      list(-500, 500),
                      data.frame(-500, 500))
    lapply(bad_input, function(x){
        expect_error(addTssAnnotation(guideSetExample,
                                      tssObject=tssObjectExample,
                                      tss_window=x))
    })
    good_input <- list(NULL,
                       c(0, 0),
                       c(-250, 50))
    lapply(good_input, function(x){
        expect_error(addTssAnnotation(guideSetExample,
                                      tssObject=tssObjectExample,
                                      tss_window=x),
                     regexp=NA)
    })
})


test_that("tss_window must be a non-positive and non-negative value pair", {
    bad_input <- list(c(-500, -100),
                      c(100, 500),
                      c(500, -500))
    lapply(bad_input, function(x){
        expect_error(addTssAnnotation(guideSetExample,
                                      tssObject=tssObjectExample,
                                      tss_window=x))
    })
    good_input <- list(c(0, 0),
                       c(-250, 50))
    lapply(good_input, function(x){
        expect_error(addTssAnnotation(guideSetExample,
                                      tssObject=tssObjectExample,
                                      tss_window=x),
                     regexp=NA)
    })
})


test_that("ignore.strand arg must be a logical value", {
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addTssAnnotation(guideSetExample[1],
                                      tssObject=tssObjectExample,
                                      ignore.strand=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(addTssAnnotation(guideSetExample,
                                      tssObject=tssObjectExample,
                                      ignore.strand=x),
                     regexp=NA)
    })
})


test_that("tssAnnotation is appended to guideSet", {
    guides <- addTssAnnotation(guideSetExample[1], tssObjectExample)
    expect_true("tssAnnotation" %in% names(mcols(guides)))
    expect_type(mcols(guides)$tssAnnotation, "S4")
    expect_true(is(mcols(guides)$tssAnnotation, "CompressedSplitDFrameList"))
})


test_that("function handles cases with no annotation to output gracefully", {
    guideSet <- guideSetExample[1]
    end(guideSet) <- start(guideSet) <- 1234
    guideSet$pam_site <- guideSet$cut_site <- 1234
    expect_error(results <- addTssAnnotation(guideSet, tssObjectExample),
                 regexp=NA)
    expect_equal(nrow(tssAnnotation(results)), 0)
})


test_that("all mcols are directly added from tssObject (excluding ID)", {
    tssMcols <- names(mcols(tssObjectExample))
    tssMcols <- setdiff(tssMcols, "ID")
    out <- addTssAnnotation(guideSetExample, tssObjectExample)
    out <- tssAnnotation(out)
    expect_true(all(tssMcols %in% colnames(out)))
})


test_that("ID column in tssObject is added as tss_id in output", {
    out <- addTssAnnotation(guideSetExample, tssObjectExample)
    out <- tssAnnotation(out)
    expect_true("tss_id" %in% colnames(out))
    expect_false("ID" %in% colnames(out))
})


test_that("appropriate anchor site is used", {
    cut_output <- addTssAnnotation(guideSetExample[1],
                                   tssObject=tssObjectExample,
                                   anchor="cut_site")
    cut_output <- tssAnnotation(cut_output)
    pam_output <- addTssAnnotation(guideSetExample[1],
                                   tssObject=tssObjectExample,
                                   anchor="pam_site")
    pam_output <- tssAnnotation(pam_output)
    expect_equal(guideSetExample$cut_site[1], cut_output$anchor_site)
    expect_equal(guideSetExample$pam_site[1], pam_output$anchor_site)
    expect_true(cut_output$anchor_site != pam_output$anchor_site)
})


test_that("annotation is only added for TSS within tss_window", {
    gs <- guideSetExample[1:2]
    
    no_annotation <- addTssAnnotation(gs, tssObjectExample, tss_window=c(0,0))
    no_annotation <- tssAnnotation(no_annotation)
    expect_equal(nrow(no_annotation), 0)
    
    one_guide_annotation <- addTssAnnotation(gs, tssObjectExample,
                                                tss_window=c(0,128))
    one_guide_annotation <- tssAnnotation(one_guide_annotation)
    expect_equal(nrow(one_guide_annotation), 1)
    
    two_guide_annotation <- addTssAnnotation(gs, tssObjectExample)
    two_guide_annotation <- tssAnnotation(two_guide_annotation)
    expect_equal(nrow(two_guide_annotation), 2)
    expect_equal(length(unique(rownames(two_guide_annotation))), 2)
    
    three_guide_annotation <- addTssAnnotation(gs, tssObjectExample,
                                               tss_window=c(-10480,500))
    three_guide_annotation <- tssAnnotation(three_guide_annotation)
    expect_equal(nrow(three_guide_annotation), 3)
    expect_equal(length(unique(rownames(three_guide_annotation))), 2)
    
    four_guide_annotation <- addTssAnnotation(gs, tssObjectExample,
                                               tss_window=c(-10500,500))
    four_guide_annotation <- tssAnnotation(four_guide_annotation)
    expect_equal(nrow(four_guide_annotation), 4)
    expect_equal(length(unique(rownames(four_guide_annotation))), 2)
})


test_that("ignore.strand argument is properly applied", {
    ignoreOn <- addTssAnnotation(guideSetExample,
                                 tssObject=tssObjectExample,
                                 ignore.strand=TRUE)
    ignoreOn <- tssAnnotation(ignoreOn)
    ignoreOff <- addTssAnnotation(guideSetExample,
                                  tssObject=tssObjectExample,
                                  ignore.strand=FALSE)
    ignoreOff <- tssAnnotation(ignoreOff)
    expect_equal(length(unique(ignoreOn$strand)), 2)
    expect_equal(length(unique(ignoreOff$strand)), 1)
    sameStrand <- ignoreOn$strand == unique(ignoreOff$strand)
    expect_equal(ignoreOn[sameStrand, , drop=FALSE], ignoreOff)
})


test_that("function does not throw error when guideSet uses custom genome", {
    guides <- findSpacers("CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTCGGGTGTTA")
    expect_error(addTssAnnotation(guides, tssObjectExample),
                 regexp=NA)
})
