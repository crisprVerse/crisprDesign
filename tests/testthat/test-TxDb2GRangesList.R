## trivial tests for getTxDb wrapper
# test_that("file argument is required to be NA or path to a GFF file", {
# })
# test_that("organism argument is required to be recognized by GenomicFeatures", {
# })
# test_that("tx_attrib argument is a valid argument for makeTxDbFromEnsembl", {
# })
# test_that("other supplied arguments are usable in makeTxDbFromEnsembl or makeTxDbFromGFF", {
# })
# test_that("getTxDb returns a TxDb object", {
# })
# test_that("getTxDb matches arguments in metadata", {
# })




test_that("txdb argument is required to be a TxDb object", {
    skip("long run time")
    
    txdb <- getTxDb()
    bad_input <- list(NULL,
                      "txdb",
                      grListExample)
    lapply(bad_input, function(x){
        expect_error(TxDb2GRangesList(x))
    })
    good_input <- list(txdb)
    lapply(good_input, function(x){
        expect_error(TxDb2GRangesList(x), regexp=NA)
    })
})


test_that("standardChromOnly argument is required to be a logical value", {
    skip("long run time")
    
    txdb <- getTxDb()
    bad_input <- list(NULL,
                      0,
                      NA,
                      "TRUE")
    lapply(bad_input, function(x){
        expect_error(TxDb2GRangesList(txdb, standardChromOnly=x))
    })
    good_input <- list(TRUE,
                       FALSE)
    lapply(good_input, function(x){
        expect_error(TxDb2GRangesList(txdb, standardChromOnly=x),
                     regexp=NA)
    })
})


test_that("genome argument is required to be a character string or NULL", {
    skip("long run time")
    
    txdb <- getTxDb()
    bad_input <- list(NA,
                      1,
                      TRUE,
                      list("hg38"),
                      c("hg38", "mm10"))
    lapply(bad_input, function(x){
        expect_error(TxDb2GRangesList(txdb, genome=x))
    })
    good_input <- list(NULL,
                       "",
                       "hg38")
    lapply(good_input, function(x){
        expect_error(TxDb2GRangesList(txdb, genome=x), regexp=NA)
    })
})


test_that("TxDb2GRangesList returns structured GRangesList object", {
    skip("long run time")
    
    txdb <- getTxDb()
    out <- TxDb2GRangesList(txdb)
    expect_true(is(out, "GRangesList"))
    regions <- c("transcripts",
                 "exons",
                 "cds",
                 "fiveUTRs",
                 "threeUTRs",
                 "introns",
                 "tss")
    expect_true(all(regions %in% names(out)))
})
