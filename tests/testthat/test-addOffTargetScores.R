test_that("guideSet argument is required to be a GuideSet object", {
    gs_as_gr <- GRanges(seqnames=seqnames(guideSetExampleFullAnnotation),
                        ranges=IRanges(start=start(guideSetExampleFullAnnotation),
                                       width=width(guideSetExampleFullAnnotation)),
                        strand=strand(guideSetExampleFullAnnotation))
    mcols(gs_as_gr) <- mcols(guideSetExampleFullAnnotation)
    names(gs_as_gr) <- names(guideSetExampleFullAnnotation)
    metadata(gs_as_gr) <- metadata(guideSetExampleFullAnnotation)
    bad_input <- list("guideSetExampleFullAnnotation",
                      as.data.frame(guideSetExampleFullAnnotation),
                      gs_as_gr)
    lapply(bad_input, function(x){
        expect_error(addOffTargetScores(x))
    })
    expect_error(addOffTargetScores(guideSetExampleFullAnnotation), regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    # not yet implemented -- need to update functions in crisprScore
    # expect_error(addOffTargetScores(guideSetExampleFullAnnotation[0]), regexp=NA)
})


test_that("guideSet is required to have spacer alignment annotation", {
    expect_error(addOffTargetScores(guideSetExample[1]))
    expect_error(addOffTargetScores(guideSetExampleFullAnnotation[1]), regexp=NA)
})


test_that("CFD and MIT scores require guideSet using SpCas9 CrisprNuclease", {
    iqsec3 <- queryTxObject(grListExample,
                            "transcripts",
                            "gene_symbol",
                            "IQSEC3")
    genome(iqsec3) <- "hg38"
    data("AsCas12a", package="crisprBase", envir=environment())
    guides <- findSpacers(
        iqsec3,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        crisprNuclease=AsCas12a)[1]
    guides <- addSpacerAlignments(
        guides,
        bowtie_index="~/crisprIndices/bowtie/hg38/hg38")
    expect_error(addOffTargetScores(guides))
})


test_that("CFD and MIT scores require spacers having lengths 19nt or 20nt", {
    guides <- guideSetExampleFullAnnotation
    guides$protospacer <- Biostrings::substr(guides$protospacer, 1, 18)
    out <- addOffTargetScores(guides)
    expect_true(all(is.na(out$score_cfd)))
    expect_true(all(is.na(out$score_mit)))
    
    guides <- guideSetExampleFullAnnotation
    nuclease <- crisprNuclease(guides)
    spacerLength(nuclease) <- 18
    metadata(guides)$CrisprNuclease <- nuclease
    expect_error(addOffTargetScores(guides))
    spacerLength(nuclease) <- 21
    metadata(guides)$CrisprNuclease <- nuclease
    expect_error(addOffTargetScores(guides))
})


test_that("max_mm argument is required to be an integer value", {
    bad_input <- list(NULL,
                      NA,
                      "1",
                      -1,
                      0.5,
                      c(1, 2))
    lapply(bad_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation[1],
                                        max_mm=x))
    })
    good_input <- list(0,
                       1,
                       100)
    lapply(good_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation[1],
                                        max_mm=x),
                     regexp=NA)
    })
})


test_that("includeDistance argument is required to be a logical value", {
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation[1],
                                        includeDistance=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation[1],
                                        includeDistance=x),
                     regexp=NA)
    })
})


test_that("offset argument is required to be a logical value", {
    bad_input <- list(NULL,
                      NA,
                      "1",
                      -1,
                      c(1, 2))
    lapply(bad_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation,
                                        offset=x))
    })
    good_input <- list(0,
                       0.5,
                       1,
                       100)
    lapply(good_input, function(x){
        expect_error(addOffTargetScores(guideSetExampleFullAnnotation,
                                        offset=x),
                     regexp=NA)
    })
})


test_that("includeDistance argument only affects MIT score", {
    out_false <- addOffTargetScores(guideSetExampleFullAnnotation,
                                    includeDistance=FALSE)
    out_true <- addOffTargetScores(guideSetExampleFullAnnotation,
                                    includeDistance=TRUE)
    expect_identical(out_false$score_cfd, out_true$score_cfd)
    expect_false(identical(out_false$score_mit, out_true$score_mit))
})


test_that("CFD and MIT scores are appended to mcols(guideSet)", {
    iqsec3 <- queryTxObject(grListExample,
                            "transcripts",
                            "gene_symbol",
                            "IQSEC3")
    genome(iqsec3) <- "hg38"
    data("SpCas9", package="crisprBase", envir=environment())
    guides <- findSpacers(
        iqsec3,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        crisprNuclease=SpCas9)[1]
    guides <- addSpacerAlignments(
        guides,
        bowtie_index="~/crisprIndices/bowtie/hg38/hg38")
    out <- addOffTargetScores(guides)
    
    score_cols <- c("score_cfd", "score_mit")
    expect_true(all(score_cols %in% colnames(mcols(out))))
})


test_that("CFD and MIT scores are appended to alignments(guideSet)", {
    iqsec3 <- queryTxObject(grListExample,
                            "transcripts",
                            "gene_symbol",
                            "IQSEC3")
    genome(iqsec3) <- "hg38"
    data("SpCas9", package="crisprBase", envir=environment())
    guides <- findSpacers(
        iqsec3,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        crisprNuclease=SpCas9)[1]
    guides <- addSpacerAlignments(
        guides,
        bowtie_index="~/crisprIndices/bowtie/hg38/hg38")
    out <- addOffTargetScores(guides)
    
    score_cols <- c("score_cfd", "score_mit")
    expect_true(all(score_cols %in% colnames(mcols(alignments(out)))))
})


test_that("aggregated scores are correctly calculated", {
    offset <- list(0,
                   0.1,
                   1,
                   100)
    lapply(offset, function(x){
        out <- addOffTargetScores(guideSetExampleFullAnnotation[3],
                                  max_mm=3,
                                  offset=x)
        agg_cfd <- 1/(sum(alignments(out)$score_cfd, na.rm=TRUE) + x)
        agg_mit <- 1/(sum(alignments(out)$score_mit, na.rm=TRUE) + x)
        names(agg_cfd) <- spacers(out)
        names(agg_mit) <- spacers(out)
        expect_equal(mcols(out)[['score_cfd']], agg_cfd)
        expect_equal(mcols(out)[['score_mit']], agg_mit)
    })
})


test_that("Off-target scores for guides with no on-targets is NA", {
    guides <- guideSetExampleFullAnnotation[3]
    offTargetIndices <- guides$alignments[[1]]$n_mismatches > 0
    guides$alignments[[1]] <- guides$alignments[[1]][offTargetIndices]
    out <- addOffTargetScores(guides)
    expect_true(is.na(out$score_cfd))
    expect_true(is.na(out$score_mit))
})


test_that("spacers alignments are annotated with correct score values", {
    out <- addOffTargetScores(guideSetExampleFullAnnotation)
    aln <- alignments(out)
    spacers <- as.character(aln$spacer)
    protospacers <- paste0(aln$protospacer, aln$pam)
    cfd_scores <- crisprScore::getCFDScores(spacers, protospacers)
    mit_scores <- crisprScore::getMITScores(spacers, protospacers)
    
    expect_identical(alignments(out)$score_cfd, cfd_scores$score)
    expect_identical(alignments(out)$score_mit, mit_scores$score)
})
