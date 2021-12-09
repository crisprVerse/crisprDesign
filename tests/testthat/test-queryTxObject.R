## tests for queryTxObject

test_that("queryTxObject only permits queries on TxDb or GRangesList object", {
    # testing Txdb input will require a GFF3 or GTF file, or connection to Ensembl server
    # expect_error(queryTxObject([TxDb object], "transcripts", "tx_id", ""),
    #              regexp=NA)
    expect_error(queryTxObject(grListExample, "transcripts", "tx_id", ""),
                 regexp=NA)
    expect_error(queryTxObject(as.list(grListExample), "transcripts", "tx_id",
                               ""))
    expect_error(queryTxObject(unlist(grListExample), "transcripts", "tx_id",
                               ""))
})


test_that("queryTxObject enforces featureType to have a specific value", {
    grListExampleNames <- as.list(c(names(grListExample), "BAD_VALUE"))
    featureTypes <- eval(formals(queryTxObject)[["featureType"]])
    lapply(grListExampleNames, function(x){
        errorType <- NULL
        if (x %in% featureTypes){
            errorType <- NA
        }
        expect_error(queryTxObject(grListExample, x, "gene_id", ""),
                     regexp=errorType)
    })
})


test_that("queryTxObject enforces queryColumn to be a character string", {
    bad_input <- list(list("tx_id"),
                      array("tx_id"),
                      data.frame("tx_id"),
                      c("tx_id", "gene_id"))
    lapply(bad_input, function(x){
        expect_error(queryTxObject(grListExample, "cds", x, ""))
    })
    expect_error(queryTxObject(grListExample, "cds", "tx_id", ""), regexp=NA)
})


test_that("queryTxObject enforces queryValue to be an atomic vector", {
    bad_input <- list(list("ENST00000538872"),
                      array("ENST00000538872"),
                      data.frame("ENST00000538872"))
    lapply(bad_input, function(x){
        expect_error(queryTxObject(grListExample, "cds", "tx_id", x))
    })
    expect_error(queryTxObject(grListExample, "cds", "tx_id",
                               "ENST00000538872"),
                 regexp=NA)
})


test_that("queryTxObject only permits queries of a single feature type", {
    bad_input <- list(c("exons", "cds"),
                      c("exons", "exons"))
    lapply(bad_input, function(x){
        expect_error(queryTxObject(grListExample, x, "tx_id", ""))
    })
    expect_error(queryTxObject(grListExample, "exons", "tx_id", ""), regexp=NA)
})


test_that("queryTxObject only permits queries of existing columns", {
    validColumns <- names(mcols(grListExample[["cds"]]))
    lapply(validColumns, function(x){
        expect_error(queryTxObject(grListExample, "cds", x, ""),
                     regexp=NA)
    })
    expect_error(queryTxObject(grListExample, "cds", "BAD_VALUE", ""))
})


test_that("queryTxObject only permits queries on a single column", {
    bad_input <- list(c("gene_id", "gene_symbol"),
                      c("gene_id", "gene_id"))
    lapply(bad_input, function(x){
        expect_error(queryTxObject(grListExample, "cds", x, ""))
    })
    expect_error(queryTxObject(grListExample, "cds", "gene_id", ""), regexp=NA)
})


test_that("queryTxObject returns empty GRanges for queryValue(s) not found", {
    expect_equal(queryTxObject(grListExample, "cds", "gene_id",
                               "NOT_A_REAL_ID"),
                 grListExample$cds[0])
})


test_that("queryTxObject queries all values in queryValue vector", {
    exon_ids <- c("ENSE00000893355", "ENSE00000893356", "NOT_A_REAL_ID")
    all_hits <- queryTxObject(grListExample, "exons", "exon_id", exon_ids)
    names(all_hits) <- NULL
    all_hits <- all_hits[order(all_hits$exon_id, all_hits$tx_id)]
    single_hits <- c(
        queryTxObject(grListExample, "exons", "exon_id", exon_ids[1]),
        queryTxObject(grListExample, "exons", "exon_id", exon_ids[2]),
        queryTxObject(grListExample, "exons", "exon_id", exon_ids[3]))
    names(single_hits) <- NULL
    single_hits <- single_hits[order(single_hits$exon_id, single_hits$tx_id)]
    expect_equal(all_hits, single_hits)
    expect_true(all(all_hits$exon_id %in% exon_ids))
    # test that no valid values in grListExample are not returned
})


test_that("queryTxObject returns unique row results", {
    gene_id <- "ENSG00000120645"
    expect_equal(queryTxObject(grListExample, "transcripts", "gene_id",
                               gene_id),
                 queryTxObject(grListExample, "transcripts", "gene_id",
                               c(gene_id, gene_id)))
})


test_that("queryTxObject returns only and all specific hits", {
    cds_len <- 3549
    cds_starts_from_query <- queryTxObject(grListExample, "exons", "cds_len",
                                      cds_len)
    cds_starts_from_query <- sort(cds_starts_from_query$cds_start)
    cds_starts_from_subset <- grListExample$exons
    cds_starts_from_subset <- sort(cds_starts_from_subset$cds_start[
        cds_starts_from_subset$cds_len == cds_len
    ])
    expect_equal(cds_starts_from_query, cds_starts_from_subset)
})


test_that("queryTxObject permits searches of NA values", {
    na_search <- queryTxObject(grListExample, "threeUTRs", "protein_id", NA)
    names(na_search) <- NULL
    three_utrs <- grListExample$threeUTRs
    names(three_utrs) <- NULL
    expect_equal(na_search, three_utrs)
})


test_that("queryTxObject permits searches of NULL and empty values", {
    expect_equal(queryTxObject(grListExample, "transcripts", "gene_id", NULL),
                 grListExample$transcripts[0])
    expect_equal(queryTxObject(grListExample, "transcripts", "gene_id", ""),
                 grListExample$transcripts[0])
})




## tests for queryTss

test_that("tssObject is required to be a GRanges objects", {
    # tssObjectExample
    # other generated GRanges
    # bad input: grListExample, data.frame, list of single gr...
    
})


test_that("queryTss enforces queryColumn to be a character string", {
    bad_input <- list(list("tx_id"),
                      # array("tx_id"), # no error
                      data.frame("tx_id"),
                      c("tx_id", "gene_id"))
    lapply(bad_input, function(x){
        expect_error(queryTss(tssObjectExample, x, ""))
    })
    expect_error(queryTss(tssObjectExample, "tx_id", ""), regexp=NA)
})


test_that("queryTss only permits queries of existing columns", {
    validColumns <- names(mcols(tssObjectExample))
    lapply(validColumns, function(x){
        expect_error(queryTss(tssObjectExample, x, ""),
                     regexp=NA)
    })
    expect_error(queryTss(tssObjectExample, "BAD_VALUE", ""))
})


test_that("queryTss only permits queries on a single column", {
    bad_input <- list(c("gene_id", "gene_symbol"),
                      c("gene_id", "gene_id"))
    lapply(bad_input, function(x){
        expect_error(queryTss(tssObjectExample, x, ""))
    })
    expect_error(queryTss(tssObjectExample, "gene_id", ""), regexp=NA)
})


test_that("queryTss enforces queryValue to be an atomic vector", {
    bad_input <- list(#list("ENST00000538872"),  # no error
                      array("ENST00000538872"))#,
                      # data.frame("ENST00000538872")) # no error
    lapply(bad_input, function(x){
        expect_error(queryTss(tssObjectExample, "tx_id", x))
    })
    expect_error(queryTss(tssObjectExample, "tx_id", "ENST00000538872"),
                 regexp=NA)
})


test_that("tss_window argument must be a integer vector of length 2 or NULL", {
    bad_input <- list(c("-500", "500"),
                      c(-500.5, 500.5),
                      c(-500),
                      c(-500, 0, 500),
                      list(-500, 500),
                      data.frame(-500, 500))
    lapply(bad_input, function(x){
        expect_error(queryTss(tssObjectExample, "gene_id", "", tss_window=x))
    })
    good_input <- list(NULL,
                       c(0, 0),
                       c(-250, 50))
    lapply(good_input, function(x){
        expect_error(queryTss(tssObjectExample, "gene_id", "", tss_window=x),
                     regexp=NA)
    })
})


test_that("tss_window must be a non-positive and non-negative value pair", {
    bad_input <- list(c(-500, -100),
                      c(100, 500),
                      c(500, -500))
    lapply(bad_input, function(x){
        expect_error(queryTss(tssObjectExample, "gene_id", "", tss_window=x))
    })
    good_input <- list(c(0, 0),
                       c(-250, 50))
    lapply(good_input, function(x){
        expect_error(queryTss(tssObjectExample, "gene_id", "", tss_window=x),
                     regexp=NA)
    })
})


test_that("queryTss returns empty GRanges for queryValue(s) not found", {
    expect_equal(queryTss(tssObjectExample, "gene_id", "NOT_A_REAL_ID"),
                 tssObjectExample[0])
})


test_that("queryTss queries all values in queryValue vector", {
    tx_ids <- c("ENST00000538872", "ENST00000382841", "NOT_A_REAL_ID")
    all_hits <- queryTss(tssObjectExample, "tx_id", tx_ids)
    names(all_hits) <- NULL
    all_hits <- all_hits[order(all_hits$tx_id)]
    single_hits <- c(
        queryTss(tssObjectExample, "tx_id", tx_ids[1]),
        queryTss(tssObjectExample, "tx_id", tx_ids[2]),
        queryTss(tssObjectExample, "tx_id", tx_ids[3]))
    names(single_hits) <- NULL
    single_hits <- single_hits[order(single_hits$tx_id)]
    expect_equal(all_hits, single_hits)
    expect_true(all(all_hits$tx_id %in% tx_ids))
    # test that no valid values in grListExample are not returned
})


test_that("queryTss returns unique row results", {
    tx_id <- "ENST00000538872"
    expect_equal(queryTss(tssObjectExample, "tx_id", tx_id),
                 queryTss(tssObjectExample, "tx_id", c(tx_id, tx_id)))
})


test_that("queryTss returns only and all specific hits", {
    id <- "IQSEC3_P1"
    results_from_query <- queryTss(tssObjectExample, "ID", id)
    names(results_from_query) <- NULL
    results_from_subset <- tssObjectExample[tssObjectExample$ID==id]
    names(results_from_subset) <- NULL
    expect_equal(results_from_query, results_from_subset)
})


test_that("queryTss permits searches of NA values", {
    na_search <- queryTss(tssObjectExample, "gene_id", NA)
    expect_equal(na_search, tssObjectExample[0])
})


test_that("queryTss permits searches of NULL and empty values", {
    expect_equal(queryTss(tssObjectExample, "gene_id", NULL),
                 tssObjectExample[0])
    expect_equal(queryTss(tssObjectExample, "gene_id", ""),
                 tssObjectExample[0])
})


test_that("queryTss sets ranges of all results to tss_window range", {
    tss_window <- c(-250, 155)
    tss_window_width <- tss_window[2] - tss_window[1]
    results <- queryTss(tssObjectExample, "gene_id", "ENSG00000120645",
                        tss_window)
    expect_equal(tss_window_width, unique(width(results)))
})

