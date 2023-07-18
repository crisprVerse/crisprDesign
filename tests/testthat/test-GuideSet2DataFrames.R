data("guideSetExample")
data("guideSetExampleFullAnnotation")


test_that("GuideSet2DataFrames returns expected output for basic GuideSet", {
    expect_error(res1 <- GuideSet2DataFrames(guideSetExample,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_error(res2 <- GuideSet2DataFrames(guideSetExample,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res3 <- GuideSet2DataFrames(guideSetExample,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res4 <- GuideSet2DataFrames(guideSetExample,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_identical(res1, res2)
    expect_identical(res3, res4)
    expect_false(identical(res1, res3))
    expect_equal(names(res1), "primary")
    expect_equal(names(res3), "primary")
    expect_type(res1, "list")
    expect_type(res3, "list")
    expect_equal(nrow(res1[["primary"]]), length(guideSetExample))
    expect_equal(nrow(res3[["primary"]]), length(guideSetExample))
})


test_that("GuideSet2DataFrames handles fully-annotated GuideSets", {
    expect_error(res1 <- GuideSet2DataFrames(guideSetExampleFullAnnotation,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_error(res2 <- GuideSet2DataFrames(guideSetExampleFullAnnotation,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res3 <- GuideSet2DataFrames(guideSetExampleFullAnnotation,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res4 <- GuideSet2DataFrames(guideSetExampleFullAnnotation,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_false(identical(res1, res2))
    expect_false(identical(res1, res3))
    expect_false(identical(res1, res4))
    expect_false(identical(res2, res3))
    expect_false(identical(res2, res4))
    expect_false(identical(res3, res4))
    expect_setequal(names(res1), c("primary", "alignments", "geneAnnotation",
                                   "tssAnnotation", "enzymeAnnotation", "snps"))
    expect_equal(names(res2), "primary")
    expect_equal(names(res3), "primary")
    expect_setequal(names(res4), c("primary", "alignments", "geneAnnotation",
                                   "tssAnnotation", "enzymeAnnotation", "snps"))
    expect_type(res1, "list")
    expect_type(res2, "list")
    expect_type(res3, "list")
    expect_type(res4, "list")
    expect_equal(nrow(res1[["primary"]]), length(guideSetExampleFullAnnotation))
    expect_equal(nrow(res2[["primary"]]), length(guideSetExampleFullAnnotation))
    expect_equal(nrow(res3[["primary"]]), length(guideSetExampleFullAnnotation))
    expect_equal(nrow(res4[["primary"]]), length(guideSetExampleFullAnnotation))
    
    ## add tests for txTable, exonTable, editedAlleles
})



test_that("GuideSet2DataFrames handles basic GuideSets with NTCs", {
    ntcs <- c("ntc_1"=paste0(rep("A", spacerLength(guideSetExample)),
                             collapse=""))
    gs <- addNtcs(guideSetExample, ntcs)
    
    expect_error(res1 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_error(res2 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res3 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res4 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_identical(res1, res2)
    expect_identical(res3, res4)
    expect_false(identical(res1, res3))
    expect_equal(names(res1), "primary")
    expect_equal(names(res3), "primary")
    expect_type(res1, "list")
    expect_type(res3, "list")
    expect_equal(nrow(res1[["primary"]]), length(gs))
    expect_equal(nrow(res3[["primary"]]), length(gs))
    
    ntc_index <- which(res1$primary$ID == "ntc_1")
    expect_equal(res1$primary$start[ntc_index], 0)
    expect_equal(res1$primary$end[ntc_index], 0)
})



test_that("GuideSet2DataFrames handles fully-annotated GuideSets with NTCs", {
    ntcs <- c("ntc_1"=paste0(
        rep("A", spacerLength(guideSetExampleFullAnnotation)),
        collapse=""))
    gs <- addNtcs(guideSetExampleFullAnnotation, ntcs)
    allAnnotationNames <- c("primary", "alignments", "geneAnnotation",
                            "tssAnnotation", "enzymeAnnotation", "snps")
    
    expect_error(res1 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_error(res2 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=TRUE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res3 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=TRUE),
                 regexp=NA)
    expect_error(res4 <- GuideSet2DataFrames(gs,
                                             useSpacerCoordinates=FALSE,
                                             primaryOnly=FALSE),
                 regexp=NA)
    expect_false(identical(res1, res2))
    expect_false(identical(res1, res3))
    expect_false(identical(res1, res4))
    expect_false(identical(res2, res3))
    expect_false(identical(res2, res4))
    expect_false(identical(res3, res4))
    expect_setequal(names(res1), allAnnotationNames)
    expect_equal(names(res2), "primary")
    expect_equal(names(res3), "primary")
    expect_setequal(names(res4), allAnnotationNames)
    expect_type(res1, "list")
    expect_type(res2, "list")
    expect_type(res3, "list")
    expect_type(res4, "list")
    expect_equal(nrow(res1[["primary"]]), length(gs))
    expect_equal(nrow(res2[["primary"]]), length(gs))
    expect_equal(nrow(res3[["primary"]]), length(gs))
    expect_equal(nrow(res4[["primary"]]), length(gs))
    
    ## add tests for txTable, exonTable, editedAlleles
})
