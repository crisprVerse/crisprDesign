# library("crisprDesignData")
# library("BSgenome.Hsapiens.UCSC.hg38")
# data("txdb_human")
# txids <- c("ENST00000607096", "ENST00000457313", "ENST00000256078")
# txids_bad <- "NOT_A_VALID_ID"
# txObject <- txdb_human
# bsgenome <- BSgenome.Hsapiens.UCSC.hg38



# test_that("'txids' argument must be a character vectors", {
#     bad_input <- list(NULL,
#                       NA,
#                       TRUE,
#                       0,
#                       as.list(txids))
#     lapply(bad_input, function(x){
#         expect_error(getMrnaSequences(txids=x,
#                                       txObject=txObject,
#                                       bsgenome=bsgenome))
#     })
#     good_input <- list(character(0),
#                        "",
#                        txids,
#                        txids_bad)
#     lapply(good_input, function(x){
#         expect_error(getMrnaSequences(txids=x,
#                                       txObject=txObject,
#                                       bsgenome=bsgenome),
#                      regexp=NA)
#     })
# })


# test_that("'txObject' argument must be a TxDb or GRangesList object", {
#     bad_input <- list(NULL,
#                       NA,
#                       "txObject",
#                       as.list(txObject),
#                       txObject[[1]])
#     lapply(bad_input, function(x){
#         expect_error(getMrnaSequences(txids=txids,
#                                       txObject=x,
#                                       bsgenome=bsgenome))
#     })
#     good_input <- list(txObject)
#     lapply(good_input, function(x){
#         expect_error(getMrnaSequences(txids=txids,
#                                       txObject=x,
#                                       bsgenome=bsgenome),
#                      regexp=NA)
#     })
# })


# test_that("'bsgenome' argument must be a BSgenome object", {
#     bad_input <- list(NULL,
#                       NA,
#                       "BSgenome.Hsapiens.UCSC.hg38",
#                       BSgenome.Hsapiens.UCSC.hg38[["chr1"]])
#     lapply(bad_input, function(x){
#         expect_error(getMrnaSequences(txids=txids,
#                                       txObject=txObject,
#                                       bsgenome=x))
#     })
#     good_input <- list(bsgenome)
#     lapply(good_input, function(x){
#         expect_error(getMrnaSequences(txids=txids,
#                                       txObject=txObject,
#                                       bsgenome=x),
#                      regexp=NA)
#     })
# })


# test_that("function returns empty DNAStringSet if no txids found", {
#     out <- getMrnaSequences(txids_bad,
#                             txObject=txObject,
#                             bsgenome=bsgenome)
#     expect_true(is(out, "DNAStringSet"))
#     expect_true(length(out) == 0)
# })


# test_that("function ignores txids not found", {
#     out <- getMrnaSequences(c(txids_bad, txids),
#                             txObject=txObject,
#                             bsgenome=bsgenome)
#     expect_true(is(out, "DNAStringSet"))
#     expect_true(length(out) == length(txids))
# })


# test_that("output is not duplicated if txids contain duplicates", {
#     out <- getMrnaSequences(c(txids, txids),
#                             txObject=txObject,
#                             bsgenome=bsgenome)
#     expect_true(is(out, "DNAStringSet"))
#     expect_equal(length(out), length(txids))
#     expect_equal(length(names(out)), length(unique(names(out))))
# })


# # update as necessary
# test_that("correct output is obtained from bsgenome", {
#     mdata <- metadata(txdb_human)
#     release <- mdata$value[mdata$name == "Ensembl release"]
#     if (release == "104"){
#         out <- getMrnaSequences(txids[1],
#                                 txObject=txObject,
#                                 bsgenome=bsgenome)
#         out <- as.character(out)[[1]]
#         # sequence obtained from ensembl.org
#         seq <- paste0("GGATGCCCAGCTAGTTTGAATTTTAGATAAACAACGAATAATTTCGTAGC",
#                       "ATAAATATGTCCCAAGCTTAGTTTGGGACATACTTATGCTAAAAAACATT",
#                       "ATTGGTTGTTTATCTGAGATTCAGAATTAAGCATTTTA")
#         expect_equal(out, seq)
#     }
# })
