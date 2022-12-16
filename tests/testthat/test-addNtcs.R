set.seed(1000)

all_ntcs <- vapply(1:10, function(x){
    seq <- sample(c("A", "C", "G", "T"), 20, replace=TRUE)
    paste0(seq, collapse="")
}, FUN.VALUE=character(1))
names(all_ntcs) <- paste0("ntc_", 1:10)



test_that("object argument must be a GuideSet object", {
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
        expect_error(addNtcs(x, all_ntcs))
    })
    expect_error(addNtcs(guideSetExample, all_ntcs), regexp=NA)
})


test_that("function handles empty object input gracefully", {
    expect_error(addNtcs(guideSetExample[0], all_ntcs),
                 regexp=NA)
})



test_that("ntcs argument must be a character vector of DNA sequences", {
    bad_input <- list(1:10,
                      TRUE,
                      as.list(all_ntcs),
                      "NOT_A_DNA_SEQUENCE")
    lapply(bad_input, function(x){
        expect_error(addNtcs(guideSetExample, x))
    })
})



test_that("ntcs argument must have unique names wrt to object and ntcs", {
    no_names <- all_ntcs
    names(no_names) <- NULL
    duplicate_names <- all_ntcs
    names(duplicate_names)[2] <- names(duplicate_names)[1]
    spacer_names <- all_ntcs
    names(spacer_names)[1] <- names(guideSetExample)[1]
    bad_input <- list(no_names,
                      duplicate_names,
                      spacer_names)
    lapply(bad_input, function(x){
        expect_error(addNtcs(guideSetExample, x))
    })
})


## test that empty ntcs returns same GuideSet




