targetOrigin <- "customSequences"

test_that("GuideSet constructor requires non-NA seqnames value(s)", {
    crisprNuclease <- SpCas9
    bad_input <- list(NULL,
                      NA,
                      list("chr1", "chr2"))
    lapply(bad_input, function(x){
        protospacers <- rep("ACTG", length(x))
        customSequences <- protospacers
        pams <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=protospacers,
                              pams=pams,
                              seqnames=x,
                              CrisprNuclease=crisprNuclease,
                              targetOrigin=targetOrigin,
                              customSequences=customSeqs))
    })
    good_input <- list("chr1",
                       1,
                       c("chr1", "chr2"),
                       "custom")
    lapply(good_input, function(x){
        protospacers <- rep("ACTG", length(x))
        customSequences <- protospacers
        pams <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=protospacers,
                              pams=pams,
                              seqnames=x,
                              CrisprNuclease=crisprNuclease,
                              targetOrigin=targetOrigin,
                              customSequences=customSequences),
                     regexp=NA)
    })
})


test_that("GuideSet constructor requires non-NA protospacer value(s)", {
    crisprNuclease <- SpCas9
    bad_input <- list(NULL,
                      NA,
                      "",
                      list("ACTG", "ACTG"))
    lapply(bad_input, function(x){
        seqnames <- rep("chr1", length(x))
        pams <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=x,
                              customSequences=x,
                              targetOrigin=targetOrigin,
                              pams=pams,
                              seqnames=seqnames,
                              CrisprNuclease=crisprNuclease))
    })
    good_input <- list("ACTG",
                       c("ACTG", "ACTG"),
                       "actg")
    lapply(good_input, function(x){
        seqnames <- rep("chr1", length(x))
        pams <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=x,
                              customSequences=x,
                              targetOrigin=targetOrigin,
                              pams=pams,
                              seqnames=seqnames,
                              CrisprNuclease=crisprNuclease),
                     regexp=NA)
    })
})


test_that("GuideSet constructor requires non-NA pam value(s)", {
    crisprNuclease <- SpCas9
    bad_input <- list(NULL,
                      NA,
                      "",
                      list("ACTG", "ACTG"))
    lapply(bad_input, function(x){
        seqnames <- rep("chr1", length(x))
        protospacers <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=protospacers,
                              customSequences=x,
                              targetOrigin=targetOrigin,
                              pams=x,
                              seqnames=seqnames,
                              CrisprNuclease=crisprNuclease))
    })
    good_input <- list("ACTG",
                       c("ACTG", "ACTG"),
                       "actg")
    lapply(good_input, function(x){
        seqnames <- rep("chr1", length(x))
        protospacers <- rep("ACTG", length(x))
        expect_error(GuideSet(protospacers=protospacers,
                              customSequences=x,
                              targetOrigin=targetOrigin,
                              pams=x,
                              seqnames=seqnames,
                              CrisprNuclease=crisprNuclease),
                     regexp=NA)
    })
})


test_that("GuideSet constructor requires a CrisprNuclease object", {
    bad_input <- list(NULL,
                      NA,
                      "SpCas9")
    lapply(bad_input, function(x){
        expect_error(GuideSet(protospacers="ACTG",
                              pams="ACTG",
                              seqnames="chr1",
                              CrisprNuclease=x,
                              customSequences="ACTG",
                              targetOrigin=targetOrigin))
    })
    expect_error(GuideSet(protospacers="ACTG",
                          pams="ACTG",
                          seqnames="chr1",
                          CrisprNuclease=SpCas9,
                          customSequences="ACTG",
                        targetOrigin=targetOrigin),
                 regexp=NA)
})


test_that("GuideSet constructor requires strand to be missing or in +/-/*", {
    bad_input <- list("pos",
                      "positive",
                      1)
    lapply(bad_input, function(x){
        expect_error(GuideSet(protospacers="ACTG",
                              customSequences="ACTG",
                              targetOrigin=targetOrigin,
                              pams="ACTG",
                              seqnames="chr1",
                              CrisprNuclease=SpCas9,
                              strand=x))
    })
    good_input <- list(NULL,
                       "+",
                       "-",
                       "*")
    lapply(good_input, function(x){
        expect_error(GuideSet(protospacers="ACTG",
                              customSequences="ACTG",
                              targetOrigin=targetOrigin,
                              pams="ACTG",
                              seqnames="chr1",
                              CrisprNuclease=SpCas9,
                              strand=x),
                     regexp=NA)
    })
    expect_warning(GuideSet(protospacers="ACTG",
                            pams="ACTG",
                            seqnames="chr1",
                            customSequences="ACTG",
                            targetOrigin=targetOrigin,
                            CrisprNuclease=SpCas9,
                            strand=NA))
})


test_that("GuideSet constructor requires pam_site to be an integer", {
    bad_input <- list(NULL,
                      NA,
                      # 1.5,     # not yet implemented
                      1i,
                      "1")
    lapply(bad_input, function(x){
        expect_error(GuideSet(protospacers="ACTG",
                              pams="ACTG",
                              seqnames="chr1",
                              customSequences="ACTG",
                              targetOrigin=targetOrigin,
                              CrisprNuclease=SpCas9,
                              pam_site=x))
    })
    good_input <- list(0,
                       -1,
                       1,
                       1L)
    lapply(good_input, function(x){
        expect_error(GuideSet(protospacers="ACTG",
                              pams="ACTG",
                              seqnames="chr1",
                              customSequences="ACTG",
                              targetOrigin=targetOrigin,
                              CrisprNuclease=SpCas9,
                              pam_site=x),
                     regexp=NA)
    })
})


test_that("GuideSet constructor stores seqinfo", {
    # gs <- GuideSet(protospacers="ACTG",
    #                pams="ACTG",
    #                seqnames="chr1",
    #                CrisprNuclease=SpCas9,
    #                genome="test")
    # expect_equal(metadata(gs)$genome, "test")
    # expect_equal(unique(genome(gs)), "test")   # not yet implemented
})


test_that("GuideSet constructor stores seqlengths", {
    
})


test_that("GuideSet consctructor stores mcol arguments", {
    
})


# different rules
# test_that("GuideSet constructor requires relevant args to be same length", {
#     lens <- c(1, 2, 4)
#     lens <- expand.grid(protospacers=lens,
#                         pams=lens,
#                         seqnames=lens,
#                         strand=lens,
#                         pam_site=lens)
#     lapply(seq_len(nrow(lens)), function(x){
#         print(lens[x,])
#         protospacers <- rep("ACTG", lens$protospacers[x])
#         pams <- rep("A", lens$pams[x])
#         seqnames <- rep("chr1", lens$seqnames[x])
#         strand <- rep("+", lens$strand[x])
#         pam_site <- rep(1, lens$pam_site[x])
#         if (length(protospacers) == length(pams) &&
#             length(protospacers) <= length(seqnames) &&
#             (length(protospacers) == length(strand) ||
#              length(strand) == 1) &&
#             (length(protospacers) == length(pam_site) ||
#              length(pam_site) == 1)){
#             regexp <- NA
#         } else {
#             regexp <- NULL
#         }
#         expect_error(GuideSet(protospacers=protospacers,
#                               pams=pams,
#                               seqnames=seqnames,
#                               pam_site=pam_site,
#                               strand=strand,
#                               CrisprNuclease=SpCas9),
#                      regexp=regexp)
#     })
# })


test_that("object argument is required to have type 'GuideSet'", {
    gs_as_gr <- GRanges(seqnames=seqnames(guideSetExample),
                        ranges=IRanges(start=start(guideSetExample),
                                       width=width(guideSetExample)))
    mcols(gs_as_gr) <- mcols(guideSetExample)
    names(gs_as_gr) <- names(guideSetExample)
    
    accessors <- list(crisprNuclease,
                      spacers,
                      pams,
                      pamSites,
                      cutSites,
                      protospacers,
                      spacerLength,
                      prototypeSequence,
                      pamLength,
                      pamSide,
                      snps,
                      alignments,
                      onTargets,
                      offTargets,
                      geneAnnotation,
                      tssAnnotation,
                      enzymeAnnotation)
    lapply(accessors, function(fun){
        expect_error(fun("guideSetExample"))
        expect_error(fun(as.data.frame(guideSetExample)))
        expect_error(fun(gs_as_gr))
        expect_error(fun(guideSetExample), regexp=NA)
    })
})


test_that("crisprNuclease accessor returns crisprNuclease object", {
    expect_true(is(crisprNuclease(guideSetExample), "CrisprNuclease"))
})


test_that("as.character TRUE returns character vector; FALSE, DNAStringSet", {
    accessors <- list(spacers,
                      pams,
                      protospacers)
    lapply(accessors, function(fun){
        expect_true(is.vector(fun(guideSetExample, as.character=TRUE),
                              mode="character"))
        expect_true(is(fun(guideSetExample, as.character=FALSE),
                       "DNAStringSet"))
    })
})


test_that("pamSites and cutSites return numeric vector of length(GuideSet)", {
    expect_equal(length(guideSetExample), length(pamSites(guideSetExample)))
    expect_equal(length(guideSetExample), length(cutSites(guideSetExample)))
    expect_true(is.vector(pamSites(guideSetExample), mode="numeric"))
    expect_true(is.vector(cutSites(guideSetExample), mode="numeric"))
})


test_that("guideSet accessor agrees with crisprNuclease accessor value", {
    accessors <- list(spacerLength,
                      prototypeSequence,
                      pamLength,
                      pamSide)
    lapply(accessors, function(fun){
        expect_equal(fun(guideSetExample), fun(SpCas9))
    })
})


test_that("annotation accessors return NULL when GuideSet lacks annotation", {
    accessors <- list(snps,
                      alignments,
                      onTargets,
                      offTargets,
                      geneAnnotation,
                      tssAnnotation,
                      enzymeAnnotation)
    lapply(accessors, function(fun){
        expect_null(fun(guideSetExample))
    })
})


test_that("unlist is appropriately applied to annotation accessors", {
    accessors_df <- list(snps,
                         geneAnnotation,
                         tssAnnotation,
                         enzymeAnnotation)
    lapply(accessors_df, function(fun){
        expect_true(is(fun(guideSetExampleFullAnnotation, unlist=TRUE),
                       "DataFrame"))
        expect_true(is(fun(guideSetExampleFullAnnotation, unlist=FALSE),
                       "DataFrameList"))
    })
    accessors_gr <- list(alignments,
                         onTargets,
                         offTargets)
    lapply(accessors_gr, function(fun){
        expect_true(is(fun(guideSetExampleFullAnnotation, unlist=TRUE),
                       "GRanges"))
        expect_true(is(fun(guideSetExampleFullAnnotation, unlist=FALSE),
                       "GRangesList"))
    })
})


test_that("alignment accessors return NULL if columnName is not valid", {
    accessors <- list(alignments,
                      onTargets,
                      offTargets)
    lapply(accessors, function(fun){
        expect_null(fun(guideSetExampleFullAnnotation,
                        columnName="BAD_COLUMN_NAME"))
        expect_null(fun(guideSetExampleFullAnnotation,
                        columnName="n0"))
    })
})


test_that("max_mismatches is appropriately applied", {
    mismatch_tally <- offTargets(guideSetExampleFullAnnotation)$n_mismatches
    mismatch_tally <- table(mismatch_tally)
    mismatch_tally <- cumsum(mismatch_tally)
    mismatch_tally <- c("0"=0, mismatch_tally, "Inf"=max(mismatch_tally))
    lapply(seq_along(mismatch_tally), function(x){
        max_mismatches <- as.numeric(names(mismatch_tally)[x])
        mismatch_count <- mismatch_tally[[x]]
        expect_equal(length(offTargets(guideSetExampleFullAnnotation,
                                       max_mismatches=max_mismatches)),
                     mismatch_count)
    })
})


test_that("gene_id filter is appropriately applied", {
    gene_ids <- geneAnnotation(guideSetExampleFullAnnotation)$gene_id
    gene_ids <- unique(gene_ids)
    lapply(gene_ids, function(x){
        id <- geneAnnotation(guideSetExampleFullAnnotation, gene_id=x)
        expect_equal(unique(id$gene_id), x)
    })
    expect_true(nrow(geneAnnotation(guideSetExampleFullAnnotation,
                                    gene_id="BAD_ID")) == 0)
})


test_that("tx_id filter is appropriately applied", {
    tx_ids <- geneAnnotation(guideSetExampleFullAnnotation)$tx_id
    tx_ids <- unique(tx_ids)
    lapply(tx_ids, function(x){
        id <- geneAnnotation(guideSetExampleFullAnnotation, tx_id=x)
        expect_equal(unique(id$tx_id), x)
    })
    expect_true(nrow(geneAnnotation(guideSetExampleFullAnnotation,
                                    tx_id="BAD_ID")) == 0)
})


test_that("gene_symbol filter is appropriately applied", {
    gene_symbols <- geneAnnotation(guideSetExampleFullAnnotation)$gene_symbol
    gene_symbols <- unique(gene_symbols)
    lapply(gene_symbols, function(x){
        id <- geneAnnotation(guideSetExampleFullAnnotation, gene_symbol=x)
        expect_equal(unique(id$gene_symbol), x)
    })
    expect_true(nrow(geneAnnotation(guideSetExampleFullAnnotation,
                               gene_symbol="BAD_SYMBOL")) == 0)
})
