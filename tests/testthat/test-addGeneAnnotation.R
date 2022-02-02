# addGeneAnnotation <- function(guideSet,
#                               txObject=NULL,
#                               anchor=c("cut_site", "pam_site"),
#                               ignore_introns=TRUE,
#                               ignore.strand=TRUE,
#                               addPfam=FALSE,
#                               mart_dataset=NULL

data("guideSetExample", "grListExample")


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
        expect_error(addGeneAnnotation(x, txObject=grListExample))
    })
    expect_error(addGeneAnnotation(guideSetExample, txObject=grListExample),
                 regexp=NA)
})


test_that("an empty guideSet is handled gracefully", {
    # not yet implemented -- problem only when adding Pfam
})


test_that("function throws error if txObject is not TxDb or GRangesList", {
    bad_input <- list(NULL,
                      grListExample[['transcripts']])
    lapply(bad_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample, txObject=x))
    })
    expect_error(addGeneAnnotation(guideSetExample, txObject=grListExample),
                 regexp=NA)
    ## long run time
    # expect_error(addGeneAnnotation(guideSetExample, txObject=getTxDb()),
    #              regexp=NA)
})


test_that("anchor argument is required to be a specific value", {
    expect_error(addGeneAnnotation(guideSetExample, txObject=grListExample,
                                   anchor="BAD_VALUE"))
    expect_error(addGeneAnnotation(guideSetExample, txObject=grListExample,
                                   anchor="cut_site"),
                 regexp=NA)
    expect_error(addGeneAnnotation(guideSetExample, txObject=grListExample,
                                   anchor="pam_site"),
                 regexp=NA)
})


test_that("ignore_introns argument is required to be a logical value", {
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample[1],
                                       txObject=grListExample,
                                       ignore_introns=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample[1],
                                       txObject=grListExample,
                                       ignore_introns=x),
                     regexp=NA)
    })
})


test_that("ignore.strand argument is required to be a logical value", {
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample[1],
                                       txObject=grListExample,
                                       ignore.strand=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample,
                                       txObject=grListExample,
                                       ignore.strand=x),
                     regexp=NA)
    })
})


test_that("addPfam argument is required to be a logical value", {
    mart_dataset <- "hsapiens_gene_ensembl"
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(addGeneAnnotation(guideSetExample[1],
                                       txObject=grListExample,
                                       addPfam=x,
                                       mart_dataset=mart_dataset))
    })
    ## long run time
    # good_input <- list(TRUE,
    #                    FALSE)
    # lapply(good_input, function(x){
    #     expect_error(addGeneAnnotation(guideSetExample[1],
    #                                    txObject=grListExample,
    #                                    addPfam=x,
    #                                    mart_dataset=mart_dataset),
    #                  regexp=NA)
    # })
})


## long run time
test_that("mart_dataset argument must be a dataset-name character string", {
    # bad_input <- list('BAD_DATASET',
    #                   c("hsapiens_gene_ensembl", "mmusculus_gene_ensembl"))
    # lapply(bad_input, function(x){
    #     expect_error(addGeneAnnotation(guideSetExample[1],
    #                                    txObject=grListExample,
    #                                    addPfam=TRUE,
    #                                    mart_dataset=x))
    # })
    # good_input <- list("hsapiens_gene_ensembl",
    #                    list("hsapiens_gene_ensembl"))
    # lapply(good_input, function(x){
    #     expect_error(addGeneAnnotation(guideSetExample[1],
    #                                    txObject=grListExample,
    #                                    addPfam=TRUE,
    #                                    mart_dataset="hsapiens_gene_ensembl"),
    #                  regexp=NA)
    # })
})


## more requirements on txObject?


test_that("geneAnnotation is appended to guideSet", {
    guides <- addGeneAnnotation(guideSetExample[1], grListExample)
    expect_true("geneAnnotation" %in% names(mcols(guides)))
    expect_type(mcols(guides)$geneAnnotation, "S4")
    expect_true(is(mcols(guides)$geneAnnotation, "CompressedSplitDFrameList"))
})


test_that("function handles cases with no annotation to output gracefully", {
    guideSet <- guideSetExample[1]
    end(guideSet) <- start(guideSet) <- 1234
    guideSet$pam_site <- guideSet$cut_site <- 1234
    expect_error(results <- addGeneAnnotation(guideSet, grListExample),
                 regexp=NA)
    expect_equal(nrow(geneAnnotation(results)), 0)
})


test_that("appropriate anchor site is used", {
    cut_output <- addGeneAnnotation(guideSetExample[1],
                                    txObject=grListExample,
                                    anchor="cut_site")
    cut_output <- geneAnnotation(cut_output)
    pam_output <- addGeneAnnotation(guideSetExample[1],
                                    txObject=grListExample,
                                    anchor="pam_site")
    pam_output <- geneAnnotation(pam_output)
    expect_equal(guideSetExample$cut_site[1], cut_output$anchor_site)
    expect_equal(guideSetExample$pam_site[1], pam_output$anchor_site)
    expect_true(cut_output$anchor_site != pam_output$anchor_site)
})


# has error in grListExample
test_that("ignore_introns removes gene annotations occuring in introns", {
    # hasInt <- addGeneAnnotation(guideSetExample,
    #                             grListExample,
    #                             ignore_introns=FALSE)
    # hasInt <- geneAnnotation(hasInt)
    # hasInt <- hasInt[!hasInt$cut_introns,]
    # noInt <- addGeneAnnotation(guideSetExample,
    #                            grListExample,
    #                            ignore_introns=TRUE)
    # noInt <- geneAnnotation(noInt)
    # all.equal(dim(hasInt), dim(noInt))
    # all.equal(hasInt, noInt)
    # 
    # hasInt[row.names(hasInt)=='spacer_486',]
    # noInt[row.names(noInt)=='spacer_486',]
    
    
})


test_that("ignore.strand argument is properly applied", {
    ignoreOn <- addGeneAnnotation(guideSetExample,
                                  txObject=grListExample,
                                  ignore.strand=TRUE)
    ignoreOn <- geneAnnotation(ignoreOn)
    ignoreOff <- addGeneAnnotation(guideSetExample,
                                   txObject=grListExample,
                                   ignore.strand=FALSE)
    ignoreOff <- geneAnnotation(ignoreOff)
    expect_equal(length(unique(ignoreOn$strand)), 2)
    expect_equal(length(unique(ignoreOff$strand)), 1)
    sameStrand <- ignoreOn$strand == unique(ignoreOff$strand)
    expect_equal(ignoreOn[sameStrand, , drop=FALSE], ignoreOff)
})


# test that addPfam works: should only test mart_dataset if/then... also, test that biomaRt is installed...
test_that("addPfam and mart_dataset arguments are properly applied", {
    # ignore tests for biomaRt installation
    expect_error(addGeneAnnotation(guideSetExample, grListExample,
                                   addPfam=FALSE, mart_dataset="BAD_DATASET"),
                 regexp=NA)
    # test that pfam column added when addPfam=TRUE, not added otherwise
})


test_that("function throws error when guideSet uses a custom genome", {
    guides <- findSpacers("CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTCGGGTGTTA")
    expect_error(addGeneAnnotation(guides, grListExample))
})


test_that("function outputs correct values", {
    out <- addGeneAnnotation(guideSetExample, grListExample) # consider subsetting
    # out_pfam <- addGeneAnnotation(guideSetExample, txObject, addPfam=TRUE,
    #                               mart_dataset="hsapiens_gene_ensembl")
    geneAnn <- geneAnnotation(out)
    geneAnnUniqueSpacer <- geneAnn[!duplicated(rownames(geneAnn)), , drop=FALSE]
    strand <- unique(as.character(strand(grListExample$transcripts)))
    
    expect_equal(length(out), length(out$geneAnnotation))
    expect_true(length(out) <= nrow(geneAnn))

    expect_true(all(rownames(geneAnn) %in% names(out)))
    expect_equal(unique(seqnames(out)), unique(geneAnn$chr))
    expect_true(all(geneAnn$pos %in% out$cut_site))
    expect_true(all(geneAnn$strand %in% c("+", "-")))
    
    # test gene_id, tx_id, protein_id
    
    # test cut_cds, cut_fiveUTRs, cut_threeUTRs: logical...always a single TRUE (if no introns)
    
    expect_false(any(geneAnn$cut_introns))
    
    sortedGeneAnn <- geneAnn[geneAnn$tx_id==geneAnn$tx_id[1], , drop=FALSE]
    sortedGeneAnn <- sortedGeneAnn[order(sortedGeneAnn$anchor_site,
                                         decreasing=(strand=="-")), , drop=FALSE]
    expect_equal(sortedGeneAnn$percentCDS, sort(sortedGeneAnn$percentCDS))
    expect_equal(sortedGeneAnn$aminoAcidIndex, sort(sortedGeneAnn$aminoAcidIndex))
    expect_equal(sortedGeneAnn$percentTx, sort(sortedGeneAnn$percentTx))
    
    # downstreamATG
    
    geneAnnRows <- vapply(out$geneAnnotation, nrow, FUN.VALUE=numeric(1),
                          USE.NAMES=FALSE)
    expect_equal(geneAnnRows, geneAnnUniqueSpacer$nIsoforms)
    geneAnnCodingRows <- vapply(out$geneAnnotation, function(x){
        sum(x$tx_id %in% grListExample$cds$tx_id)
    }, FUN.VALUE=numeric(1), USE.NAMES=FALSE)
    expect_equal(geneAnnCodingRows, geneAnnUniqueSpacer$nCodingIsoforms)
    
    
    isoformTestCols <- c("gene_id", "totalIsoforms", "totalCodingIsoforms")
    expect_equal(length(unique(geneAnn$gene_id)),
                 nrow(unique(geneAnn[, isoformTestCols, drop=FALSE])))
    
    expect_true(all(geneAnn$nIsoforms <= geneAnn$totalIsoforms))
    expect_true(all(geneAnn$nCodingIsoforms <= geneAnn$totalCodingIsoforms))
    expect_equal(geneAnn$percentIsoforms,
                 round(geneAnn$nIsoforms/geneAnn$totalIsoforms*100,1))
    expect_equal(geneAnn$percentCodingIsoforms,
                 round(geneAnn$nCodingIsoforms/geneAnn$totalCodingIsoforms*100,
                       1))
    
    expect_equal(geneAnn$isCommonExon,
                 geneAnn$nIsoforms == geneAnn$totalIsoforms)
    expect_equal(geneAnn$isCommonCodingExon,
                 geneAnn$nCodingIsoforms == geneAnn$totalCodingIsoforms)
})
