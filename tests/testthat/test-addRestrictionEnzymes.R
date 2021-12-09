data(restrictionEnzymes, package="crisprBase")

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
        expect_error(addRestrictionEnzymes(x))
        expect_error(getRestrictionEnzymes(x))
    })
    expect_error(addRestrictionEnzymes(guideSetExample), regexp=NA)
    expect_error(getRestrictionEnzymes(guideSetExample), regexp=NA)
})


test_that("empty guideSet input is handled gracefully", {
    expect_error(addRestrictionEnzymes(guideSetExample[0],
                                       enzymeNames='PhoI',
                                       patterns=c(enz1='ARGNT')),
                 regexp=NA)
    expect_error(getRestrictionEnzymes(guideSetExample[0],
                                       enzymeNames='PhoI',
                                       patterns=c(enz1='ARGNT')),
                 regexp=NA)
})


test_that("enzymeNames argument requires a character vector or NULL", {
    bad_input <- list(0,
                      NA,
                      restrictionEnzymes[[1]],
                      list(names(restrictionEnzymes[[1]])))
    good_input <- list(NULL,
                       character(0),
                       names(restrictionEnzymes)[1],
                       "PhoI")
    lapply(bad_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1],
                                           enzymeNames=x,
                                           includeDefault=FALSE))
        expect_error(getRestrictionEnzymes(guideSetExample[1],
                                           enzymeNames=x,
                                           includeDefault=FALSE))
    })
    lapply(good_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1], enzymeNames=x),
                     regexp=NA)
        expect_error(getRestrictionEnzymes(guideSetExample[1], enzymeNames=x),
                     regexp=NA)
    })
})


test_that("enzymeNames values must be in crisprBase::restrictionEnzymes", {
    expect_error(addRestrictionEnzymes(guideSetExample[1],
                                       enzymeNames="BAD_NAME"))
    expect_error(addRestrictionEnzymes(guideSetExample[1],
                                       enzymeNames=names(
                                           restrictionEnzymes)[1]),
                 regexp=NA)
})


test_that("patterns argument requires a named character vector or NULL", {
    bad_input <- list(restrictionEnzymes[[1]],
                      c(a=1, b=2),
                      c('ACGT', 'AGCT'),
                      c(a='ACGT', 'AGCT'))
    lapply(bad_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1], patterns=x))
        expect_error(getRestrictionEnzymes(guideSetExample[1], patterns=x))
    })
    expect_error(addRestrictionEnzymes(guideSetExample[1],
                                       patterns=c(a='ACGT', b='AGCT')),
                 regexp=NA)
    expect_error(getRestrictionEnzymes(guideSetExample[1],
                                       patterns=c(a='ACGT', b='AGCT')),
                 regexp=NA)
})


test_that("patterns argument is required to have unique names", {
    # currently, function silently ignores duplicates
})


test_that("patterns argument is required to have only IUPAC characters", {
    bad_input <- list(c(enz1="ACFG"),
                      c(enz1="AC-G"))
    lapply(bad_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample, patterns=x))
        expect_error(getRestrictionEnzymes(guideSetExample, patterns=x))
    })
    good_input <- list(c(enz1="ACNG"))
    lapply(good_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample, patterns=x),
                     regexp=NA)
        expect_error(getRestrictionEnzymes(guideSetExample, patterns=x),
                     regexp=NA)
    })
})


test_that("includeDefault argument requires a logical value", {
    bad_input <- list(NA,
                      "T")
    lapply(bad_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1],
                                           includeDefault=x,
                                           enzymeNames="EcoRI"))
        expect_error(getRestrictionEnzymes(guideSetExample[1],
                                           includeDefault=x,
                                           enzymeNames="EcoRI"))
    })
    good_input <- list(TRUE,
                       FALSE,
                       "TRUE",
                       0)
    lapply(good_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1],
                                           includeDefault=x,
                                           enzymeNames="EcoRI"),
                     regexp=NA)
        expect_error(getRestrictionEnzymes(guideSetExample[1],
                                           includeDefault=x,
                                           enzymeNames="EcoRI"),
                     regexp=NA)
    })
})


test_that("flanking5 and flanking3 arguments require DNA (character) string", {
    bad_input <- list("QWERTY",
                      "ACGN",
                      c("ACGT", "ACGT"))
    good_input <- list(NULL,
                       "",
                       "ACGT",
                       "ACGU",
                       "acgt")
    lapply(bad_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1], flanking5=x))
        expect_error(addRestrictionEnzymes(guideSetExample[1], flanking3=x))
        expect_error(getRestrictionEnzymes(guideSetExample[1], flanking5=x))
        expect_error(getRestrictionEnzymes(guideSetExample[1], flanking3=x))
    })
    lapply(good_input, function(x){
        expect_error(addRestrictionEnzymes(guideSetExample[1], flanking5=x),
                     regexp=NA)
        expect_error(addRestrictionEnzymes(guideSetExample[1], flanking3=x),
                     regexp=NA)
        expect_error(getRestrictionEnzymes(guideSetExample[1], flanking5=x),
                     regexp=NA)
        expect_error(getRestrictionEnzymes(guideSetExample[1], flanking3=x),
                     regexp=NA)
    })
})


test_that("at least one restriction enzyme input is required", {
    expect_error(addRestrictionEnzymes(guideSetExample, includeDefault=FALSE))
    expect_error(getRestrictionEnzymes(guideSetExample, includeDefault=FALSE))
})


test_that("includeDefault adds defaults restriction enzymes to output", {
    enzymes <- c("EcoRI", "KpnI", "BsmBI", "BsaI", "BbsI", "PacI", "MluI")
    guides <- addRestrictionEnzymes(guideSetExample[1])
    expect_true(all(enzymes %in% colnames(enzymeAnnotation(guides))))
})


test_that("spacers matching sites for includeDefault are correctly identified", {
    # guide for each positive enzyme case
})


test_that("spacers matching sites for enzymeNames are correctly identified", {
    # test with normal ACGT bases
    # test with non-ACGT symbols
})


test_that("spacers matching sites for patterns are correctly identified", {
    # test with normal ACGT bases
    # test with non-ACGT symbols
    # test with NNNN (should match all)
})


test_that("restriction site recognition incorporates flanking sequences", {
    # test that uses default flanking (both) (overlap both flanking and seq)
    # test using custom flanking (overlap both flanking and seq)
})
