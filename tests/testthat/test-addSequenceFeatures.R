data(SpCas9, package="crisprBase")

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
        expect_error(addSequenceFeatures(x))
    })
    expect_error(addSequenceFeatures(guideSetExample), regexp=NA)
})


test_that("function handles empty GuideSet input gracefully", {
    expect_error(addSequenceFeatures(guideSetExample[0]),
                 regexp=NA)
})


test_that("addHairpin argument must be a logical value", {
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addSequenceFeatures(guideSetExample[1],
                                         addHairpin=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(addSequenceFeatures(guideSetExample[1],
                                         addHairpin=x),
                     regexp=NA)
    })
})


test_that("backbone argument must be a nucleotide character string", {
    bad_input <- list(NULL,
                      NA,
                      "ASDF",
                      c("ACGT", "GACT"),
                      guideSetExample)
    lapply(bad_input, function(x){
        expect_error(addSequenceFeatures(guideSetExample[1], addHairpin=TRUE,
                                         backbone=x))
    })
    good_input <- list("",
                       "AGGCTAGTCCGT",
                       list("AGGCTAGTCCGT"),
                       "AGGCUAGUCCGU")
    lapply(good_input, function(x){
        expect_error(addSequenceFeatures(guideSetExample[1], addHairpin=TRUE,
                                         backbone=x),
                     regexp=NA)
    })
})


test_that("sequence feature columns are appended to GuideSet object", {
    feature_cols <- c("polyA", "polyC", "polyG", "polyT", "startingGGGGG")
    out <- addSequenceFeatures(guideSetExample[1])
    expect_true(all(feature_cols %in% colnames(mcols(out))))
})


test_that("hairpin columns only added when addHairpin=TRUE", {
    hairpin_cols <- c("selfHairpin", "backboneHairpin")
    no_hairpin <- addSequenceFeatures(guideSetExample[1], addHairpin=FALSE)
    expect_false(any(hairpin_cols %in% colnames(mcols(no_hairpin))))
    with_hairpin <- addSequenceFeatures(guideSetExample[1], addHairpin=TRUE)
    expect_true(all(hairpin_cols %in% colnames(mcols(with_hairpin))))
})


test_that("percentGC is correctly calculated", {
    spacer_input <- c("ATATATATATATATATATAT", "GTATATATATATATATATAT",
                      "GCATATATATATATATATAT", "GCGTATATATATATATATAT",
                      "GCGCATATATATATATATAT", "GCGCGTATATATATATATAT",
                      "GCGCGCATATATATATATAT", "GCGCGCGTATATATATATAT",
                      "GCGCGCGCATATATATATAT", "GCGCGCGCGTATATATATAT",
                      "GCGCGCGCGCATATATATAT", "GCGCGCGCGCGTATATATAT",
                      "GCGCGCGCGCGCATATATAT", "GCGCGCGCGCGCGTATATAT",
                      "GCGCGCGCGCGCGCATATAT", "GCGCGCGCGCGCGCGTATAT",
                      "GCGCGCGCGCGCGCGCATAT", "GCGCGCGCGCGCGCGCGTAT",
                      "GCGCGCGCGCGCGCGCGCAT", "GCGCGCGCGCGCGCGCGCGT",
                      "GCGCGCGCGCGCGCGCGCGC")
    guides <- GuideSet(protospacers=spacer_input,
                       customSequences=spacer_input,
                       targetOrigin="customSequences",
                       pams=rep("CGG", 21),
                       seqnames=rep("chr1", 21),
                       CrisprNuclease=SpCas9,
                       strand="+")
    out <- addSequenceFeatures(guides)
    expect_equal(out$percentGC, seq(0, 100, by=5))
})


test_that("spacers with homopolymers are correctly identified", {
    spacer_input <- c("CGGAAAACTTGCATGATATG",
                      "GAGTCAGCTTGCCCCATATG",
                      "CATGCAGCTTGGGGTATATG",
                      "CATGCATTTTGACGTCTATG")
    guides <- GuideSet(ids=seq_along(spacer_input),
                       protospacers=spacer_input,
                       customSequences=spacer_input,
                       targetOrigin="customSequences",
                       pams=rep("CGG", 4),
                       seqnames=rep("chr1", 4),
                       CrisprNuclease=SpCas9,
                       strand="+")
    out <- addSequenceFeatures(guides)
    expect_equal(which(out$polyA), 1)
    expect_equal(which(out$polyC), 2)
    expect_equal(which(out$polyG), 3)
    expect_equal(which(out$polyT), 4)
})


test_that("spacers with startingGGGGG are correctly identified", {
    spacer_input <- c("GGGGGAGCTTGCATGATATG",
                      "GGGGCAGCTTGCATGATATG",
                      "CATGCAGCTTGGGGGATATG")
    hasStartingGGGGG <- c(TRUE, FALSE, FALSE)
    guides <- GuideSet(ids=seq_along(spacer_input),
                       protospacers=spacer_input,
                       customSequences=spacer_input,
                       targetOrigin="customSequences",
                       pams=rep("CGG", 3),
                       seqnames=rep("chr1", 3),
                       CrisprNuclease=SpCas9,
                       strand="+")
    startingGGGGG <- addSequenceFeatures(guides)$startingGGGGG
    expect_equal(startingGGGGG, hasStartingGGGGG)
})


test_that("selfHairpin and backboneHairpin are correctly predicted", {
    spacers <- c("GAATccaagcgggacATTCa"=FALSE,    # percentGC < 50
                 "GACTccaagcgggacAGTCa"=TRUE,     # percentGC >= 50
                 "GACTccaAGTCaagcgggac"=FALSE,    # loop length 3
                 "GACTccaaAGTCagcgggac"=TRUE,     # loop length 4
                 "GACaccaaaGTCagcgggac"=FALSE,    # stem length 3
                 "GACAccaaTGTCagcgggac"=TRUE,     # stem length 4
                 "ACAccaaTGTCagcgggaca"=TRUE,     # with flanking 5' G
                 "ACATccaacaggagcgggAT"=TRUE)     # with flanking 3' GT
    guides <- GuideSet(ids=seq_along(spacers),
                       protospacers=names(spacers),
                       customSequences=names(spacers),
                       targetOrigin="customSequences",
                       pams=rep("CGG", 8),
                       seqnames=rep("chr", 8),
                       CrisprNuclease=SpCas9,
                       strand="+")
    out <- addSequenceFeatures(guides, addHairpin=TRUE)
    expect_equal(mcols(out)$selfHairpin, as.logical(spacers))
    
    spacer <- "gCTAGTGGagaggaggaggg"
    backbones <- list("aggACTAGggaa"=FALSE,       # percent GC < 50
                      "agCCACTgagga"=TRUE,        # percent GC >= 50
                      "aggACTAaggaa"=FALSE,       # stem length 4
                      "agCCACTgagga"=TRUE,        # stem length 5
                      "agaTAGCCggaa"=TRUE,        # with flanking 5' G
                      "aggACCCCaggaa"=TRUE)       # with flanking 3' GT
    guides <- GuideSet(ids=seq_along(spacer),
                       protospacers=spacer,
                       customSequences=spacer,
                       targetOrigin="customSequences",
                       pams="CGG",
                       seqnames="chr",
                       CrisprNuclease=SpCas9,
                       strand="+")
    lapply(seq_along(backbones), function(x){
        out <- addSequenceFeatures(guides, addHairpin=TRUE,
                                   backbone=names(backbones)[x])
        expect_equal(as.logical(out$backboneHairpin), backbones[[x]])
    })
})
