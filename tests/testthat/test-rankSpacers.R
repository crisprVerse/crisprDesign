update <- TRUE

gs_as_gr <- GRanges(seqnames=seqnames(guideSetExample),
                    ranges=IRanges(start=start(guideSetExample),
                                   width=width(guideSetExample)))
mcols(gs_as_gr) <- mcols(guideSetExample)
names(gs_as_gr) <- names(guideSetExample)

gs_as_gr_full <- GRanges(seqnames=seqnames(guideSetExampleFullAnnotation),
                    ranges=IRanges(start=start(guideSetExampleFullAnnotation),
                                   width=width(guideSetExampleFullAnnotation)))
mcols(gs_as_gr_full) <- mcols(guideSetExampleFullAnnotation)
names(gs_as_gr_full) <- names(guideSetExampleFullAnnotation)






test_that("validCriteria requires GuideSet input", {
    errorMessage <- "Object must be a GuideSet" # temporary check
    expect_error(validCriteria("guideSetExample"),
                 regexp=errorMessage)
    expect_error(validCriteria(as.data.frame(guideSetExample)),
                 regexp=errorMessage)
    expect_error(validCriteria(gs_as_gr),
                 regexp=errorMessage)
    expect_error(validCriteria(guideSetExampleFullAnnotation),
                 regexp=NA)
})



test_that('validCriteria gives correct output', {
    # expect_equal_to_reference(validCriteria(guideSetExample),
    #                           update=update,
    #                           file=file.path("objects/results_validcriteria_gs_no_annotation.rds"))
    # expect_equal_to_reference(validCriteria(guideSetExampleFullAnnotation),
    #                           update=update,
    #                           file=file.path("objects/results_validcriteria_gs_full_annotation_2mm.rds"))
})



test_that('validCriteria output is independent of GuideSet length', {
    expect_equal(validCriteria(guideSetExample),
                 validCriteria(guideSetExample[0]))
    expect_equal(validCriteria(guideSetExampleFullAnnotation),
                 validCriteria(guideSetExampleFullAnnotation[0]))
})





test_that("filterSpacers and rankSpacers require GuideSet input", {
    criteria <- list(polyT=FALSE)
    # errorMessage <- "Object must be a GuideSet"
    guideSets <- list("guideSetExampleFullAnnotation",
                      as.data.frame(guideSetExampleFullAnnotation),
                      gs_as_gr_full)
    lapply(guideSets, function(x){
        expect_error(filterSpacers(x, criteria), regexp=NULL)
        expect_error(rankSpacers(x, criteria), regexp=NULL)
    })
    expect_error(filterSpacers(guideSetExampleFullAnnotation, criteria),
                 regexp=NA)
    expect_error(rankSpacers(guideSetExampleFullAnnotation, criteria),
                 regexp=NA)
})



test_that("criteria is required to be a named list with non-empty names", {
    # errorMessage <- "criteria must be a non-empty named list"
    criteria_with_missing_name <- list(FALSE, 1)
    names(criteria_with_missing_name) <- "polyT"
    bad_criteria <- list("polyT=FALSE",
                         list(),
                         list(FALSE, 1),
                         list(polyT=FALSE, n1=1)[0],
                         list(polyT=FALSE, 1),
                         criteria_with_missing_name)
    good_criteria <- list(data.frame(polyT=FALSE, n1=1),
                          list(polyT=FALSE, n1=1))
    lapply(bad_criteria, function(x){
        expect_error(filterSpacers(guideSetExampleFullAnnotation, x, regexp=NULL))
        expect_error(rankSpacers(guideSetExampleFullAnnotation, x, regexp=NULL))
    })
    lapply(good_criteria, function(x){
        expect_error(filterSpacers(guideSetExampleFullAnnotation, x, regexp=NA))
        expect_error(rankSpacers(guideSetExampleFullAnnotation, x, regexp=NA))
    })
})




criteria_no_id <- list(polyT=TRUE)
criteria_txId_only <- list(cut_cds=TRUE)
criteria_geneId_only <- list(isCommonExon=TRUE)
criteria_either_id <- list(percentCDS=82)      # also test that all isoformAggFun work as intended

txId_bad_type <- list("ENST00000538872")
txId_bad_value <- "not_a_valid_tx_id"
txId_bad_length <- c("ENST00000538872", "ENST00000382841")
txId_valid <- "ENST00000538872"
geneId_bad_type <- list("ENSG00000120645")
geneId_bad_value <- "not_a_valid_gene_id"
geneId_bad_length <- c("ENSG00000120645", "ENSG00000256694")
geneId_valid <- "ENSG00000120645"



test_that("neither txId nor geneID are validated if not used", {
    criteria <- list(polyT=TRUE)
    txIds <- list(NULL,
                  txId_bad_type,
                  txId_bad_value,
                  txId_bad_length)
    geneIds <- list(NULL,
                    geneId_bad_type,
                    geneId_bad_value,
                    geneId_bad_length)
    lapply(seq_along(txIds), function(x){
        expect_error(filterSpacers(guideSetExampleFullAnnotation, criteria,
                                  txId=txIds[x], geneId=geneIds[x]),
                     regexp=NA)
        expect_error(rankSpacers(guideSetExampleFullAnnotation, criteria,
                                txId=txIds[x], geneId=geneIds[x]),
                     regexp=NA)
    })
})



test_that("txId is validated when used", {
    # error_bad_type <- "txId must be a length 1 character vector"
    # error_bad_value <- "txId not found in guideSet:"
    # error_gene_id <- "Supplied criteria requires a geneId"
    criterias <- list(criteria_txId_only, criteria_either_id) # for either_id/both, include cases where both IDs are provided
    # add geneId_only?
    geneIds <- list(NULL,
                    txId_bad_type,
                    txId_bad_value,
                    txId_bad_length) # add txId_valid that gives no error
    for (i in criterias){
        for (ii in geneIds){
            expect_error(filterSpacers(guideSetExampleFullAnnotation, criteria=i, txId=ii),
                         regexp=NULL)
            expect_error(rankSpacers(guideSetExampleFullAnnotation, criteria=i, txId=ii),
                         regexp=NULL)
        }
    }
})



test_that("txId is used to obtain GuideSet values correctly", {
    # expect_equal_to_reference(
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_txId_only, txId=txId_valid),
    #     update=update,
    #     file=file.path("objects/results_rankSpacers_txId_only.rds"))
    # expect_equal_to_reference(
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, txId=txId_valid),
    #     update=update,
    #     file=file.path("objects/results_rankSpacers_either_id_txId_input.rds"))
    expect_equal(
        rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, txId=txId_valid),
        rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, txId=txId_valid, geneId=geneId_valid))
})



test_that("geneId is validated when used", {
    criterias <- list(criteria_geneId_only, criteria_either_id) # for either_id/both, include cases where both IDs are provided
    # add txId_only?
    geneIds <- list(NULL,
                    geneId_bad_type,
                    geneId_bad_value,
                    geneId_bad_length) # add geneId_valid that gives no error
    for (i in criterias){
        for (ii in geneIds){
            expect_error(filterSpacers(guideSetExampleFullAnnotation, criteria=i, geneId=ii),
                         regexp=NULL)
            expect_error(rankSpacers(guideSetExampleFullAnnotation, criteria=i, geneId=ii),
                         regexp=NULL)
        }
    }
})



test_that("geneId is used to obtain GuideSet values correctly", {
    # expect_equal_to_reference(
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_geneId_only, geneId=geneId_valid),
    #     update=update,
    #     file=file.path("objects/results_rankSpacers_geneId_only.rds"))
    # expect_equal_to_reference(
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, geneId=geneId_valid),
    #     update=update,
    #     file=file.path("objects/results_rankSpacers_either_id_geneId_input.rds"))
    # expect_false(identical(
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, geneId=geneId_valid),
    #     rankSpacers(guideSetExampleFullAnnotation, criteria_either_id, txId=txId_valid, geneId=geneId_valid)))
})


test_that("isoformAggFun aggregates gene-level values as intended", {
    # tests for: max, mean, min, other (raises error)
})



test_that("for criteria attributes not found in guideSet, user is informed how to add missing annotation", {
    # test 1 column for each annotation type in gs_no_annotation
    # geneAnnotation, offTargetScores, onTargetScores, pamScores, repeats, enzymeAnnotation,
    #   sequenceFeatures, snpAnnotation, spacerAlignments, tssAnnotation
    # test for criteria not found anywhere in guideSet, ex: "nonexistent_spacer_attribute"
})



test_that("values to bin are obtainable for all guideSet annotation types", {
    # geneAnnotation, offTargetScores, onTargetScores, pamScores, repeats, enzymeAnnotation,
    #   sequenceFeatures, snpAnnotation, spacerAlignments, tssAnnotation
    # for geneAnnotation, test both txId and geneId (values that have NAs)
})


test_that("constraints for logical criteria values are enforced", {
    # criteria with length!=1, !is.logical, and !is.vector
})

test_that("logical criteria values are appropriately binned", {
    # mcol logical, geneAnnotation logical (with NAs)
    # criteria with duplicates, ex: list(polyT=FALSE, polyT=TRUE), and geneAnnotation
})


test_that("constraints for asc criteria values are enforced", {
    # criteria that are !is.numeric, !is.vector, in bad order, has duplicates
})


test_that("asc criteria values are appropriately binned", {
    # criteria with -Inf, Inf, both
    # guideSet with NA values (geneAnnotation, aln summary)
    # criteria with single bin value
    # criteriawith all bins outside of range
    # normal, detailed binning, ex: 1, 2, 3, ... (also tests edge cases)
})


test_that("constraints for desc criteria values are enforced", {
    # criteria that are !is.numeric, !is.vector, in bad order, has duplicates
})


test_that("desc criteria values are appropriately binned", {
    # criteria with -Inf, Inf, both
    # guideSet with NA values (geneAnnotation, score?)
    # criteria with single bin value
    # criteriawith all bins outside of range (ex: negative values)
    # normal, detailed binning, ex: 0.9, 0.8, 0.7, ... (also tests edge cases)
})


test_that("constraints for ranged criteria values are enforced", {
    # criteria that is !numeric, length!=2, outside of 0-100; 0-100 exactly, desc, duplicated, normal
})



test_that("ranged criteria values are appropriately binned", {
    # criteria with duplicated values with: no match, some match
    # full range; normal range
})



test_that("criteria containing duplicate attributes are handled appropriately", {
    # criteria with 2- and 3-repeats (check for ranking order and colnames)
})



test_that("output is not dependent on guideSet names", {
    # input guideSet that is unordered; has unusual/custom names; has duplicated names?
})



test_that("rankSpacers adds ranks appropriately and orders GuideSet by rank", {
    # criteria with single element...; criteria with all data types (and both logical variations)
    # test to compare ranked output against itself: 
    #       x <- rankSpacers(...), expect_equal(x, rankSpacers(x)))
})


test_that("filterSpacers only retains guides having rank==1", {
    # normal criteria with expected guides passing
    # criteria with no guides passing
    # criteria where no guides can pass (ex: list(polyT=TRUE, polyT=FALSE))
})


test_that("empty GuideSet input is returned unchanged", {
    criteria <- list(polyT=FALSE, n1=c(0, 1, 5)) # should include one from each annotation type
    expect_equal(guideSetExampleFullAnnotation[0],
                 filterSpacers(guideSetExampleFullAnnotation[0], criteria))
    expect_equal(guideSetExampleFullAnnotation[0],
                 rankSpacers(guideSetExampleFullAnnotation[0], criteria))
})
