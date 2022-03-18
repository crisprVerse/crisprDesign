# library(crisprDesignData)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Mmusculus.UCSC.mm10)
# data("txdb_human")
# data("guideSetExample", "grListExample")
# data("SpCas9", "AsCas12a")
# gs <- head(guideSetExample)
# spacers <- protospacers(gs, as.character=TRUE)


# bsgenome_human <- BSgenome.Hsapiens.UCSC.hg38
# bsgenome_mouse <- BSgenome.Mmusculus.UCSC.mm10
# #bowtie_index <- "~/crisprIndices/bowtie/hg38/hg38"
# #bwa_index <- "~/crisprIndices/bwa/hg38/hg38"
# #bowtie_index_mouse <- "~/crisprIndices/bowtie/mm10/mm10"

# custom_seq1 <- getSeq(bsgenome_human, names="chr12", start=66000, end=70000)
# custom_seq2 <- "ACTG"
# custom_seq3 <- paste0(protospacers(gs), pams(gs), collapse="")
# custom_seq4 <- .revCompBs(custom_seq3)
# custom_seq5 <- DNAStringSet(list(custom_seq1, bsgenome_human[["chrM"]]))



# test_that("getSpacerAlignments requires GuideSet or spacer input", {
#     # bad_input <- list(NULL,
#     #                   NA,
#     #                   0,
#     #                   TRUE,
#     #                   as.data.frame(gs),
#     #                   as(gs, "GRanges"))
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=x,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      crisprNuclease=SpCas9))
#     # })
#     ## skipped due to long runtime
#     # good_input <- list(gs,
#     #                    spacers,
#     #                    DNAStringSet(spacers),
#     #                    DNAString(spacers[1]))
#     # lapply(good_input, function(x){
#     #     expect_error(
#     #         suppressWarnings(
#     #             getSpacerAlignments(spacers=x,
#     #                                 aligner="bowtie",
#     #                                 aligner_index=bowtie_index,
#     #                                 bsgenome=bsgenome_human,
#     #                                 crisprNuclease=SpCas9)
#     #             ),
#     #         regexp=NA)
#     # })
# })


# test_that("getSpacerAlignments extracts spacers from GuideSet object", {
#     ## skipped due to long runtime
#     # out_gs <- getSpacerAlignments(spacers=gs,
#     #                               aligner="bowtie",
#     #                               aligner_index=bowtie_index,
#     #                               bsgenome=bsgenome_human,
#     #                               crisprNuclease=SpCas9)
#     # out_spacers <- getSpacerAlignments(spacers=spacers,
#     #                                    aligner="bowtie",
#     #                                    aligner_index=bowtie_index,
#     #                                    bsgenome=bsgenome_human,
#     #                                    crisprNuclease=SpCas9)
#     # expect_equal(out_gs, out_spacers)
# })


# test_that("aligner argument must be one of given options", {
#     # bad_input <- list(NA,
#     #                   "",
#     #                   0)
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner=x,
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      crisprNuclease=SpCas9))
#     # })
# })


# test_that("custom_seq must be valid character vector or XString[Set] object", {
#     # bad_input <- list(NULL,
#     #                   NA,
#     #                   0,
#     #                   bsgenome_human,
#     #                   list(custom_seq1, custom_seq2),
#     #                   "BAD_SEQUENCE")
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="biostrings",
#     #                                      crisprNuclease=SpCas9,
#     #                                      custom_seq=x))
#     # })
#     ## skipped due to long runtime
#     # good_input <- list(custom_seq1,
#     #                    custom_seq2)
#     # lapply(good_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="biostrings",
#     #                                      crisprNuclease=SpCas9,
#     #                                      custom_seq=x),
#     #                  regexp=NA)
#     # })
# })


# test_that("custom_seq names are appropriately assigend to found alignments", {
#     named_seq <- custom_seq5
#     names(named_seq) <- c("test1", "test2")
#     aln <- getSpacerAlignments(spacers=spacers,
#                                aligner="biostrings",
#                                crisprNuclease=SpCas9,
#                                custom_seq=named_seq)
#     aln_seqlevels <- as.character(seqlevels(aln))
#     expect_equal(sort(aln_seqlevels), sort(names(named_seq)))
# })



# test_that("aligner_index must be appropriate for selected aligner", {
#     # expect_error(getSpacerAlignments(spacers=spacers,
#     #                                  aligner="bowtie",
#     #                                  aligner_index=bwa_index,
#     #                                  bsgenome=bsgenome_human,
#     #                                  crisprNuclease=SpCas9))
#     # expect_error(getSpacerAlignments(spacers=spacers,
#     #                                  aligner="bwa",
#     #                                  aligner_index=bowtie_index,
#     #                                  bsgenome=bsgenome_human,
#     #                                  crisprNuclease=SpCas9))
# })


# test_that("aligner_index files must be appropriate for provided BSgenome", {
#     # expect_error(getSpacerAlignments(spacers=spacers,
#     #                                  aligner="bowtie",
#     #                                  aligner_index=bowtie_index,
#     #                                  bsgenome=bsgenome_mouse,
#     #                                  crisprNuclease=SpCas9))
# })


# test_that("bsgenome arg must be a BSgenome object for indexed alignments", {
#     # bad_input <- list(NULL,
#     #                   NA,
#     #                   "bsgenome_human",
#     #                   bsgenome_human[[1]])
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      crisprNuclease=SpCas9,
#     #                                      bsgenome=x))
#     # })
# })


# test_that("n_mismatches arg must be a single non-negative integer value", {
#     # bad_input <- list(NA,
#     #                   -1,
#     #                   0.5,
#     #                   "1",
#     #                   list(1),
#     #                   c(1,1))
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      crisprNuclease=SpCas9,
#     #                                      bsgenome=bsgenome_human,
#     #                                      n_mismatches=x))
#     # })
# })


# test_that("n_mismatches limit for bowtie is enforced", {
#     # bowtie_limit <- 3
#     # expect_error(getSpacerAlignments(spacers=spacers,
#     #                                  aligner="bowtie",
#     #                                  aligner_index=bowtie_index,
#     #                                  crisprNuclease=SpCas9,
#     #                                  bsgenome=bsgenome_human,
#     #                                  n_mismatches=bowtie_limit+1))
# })


# test_that("n_mismatches arg is appropriately applied", {
#     ## skip due to long runtime
#     # lapply(0:3, function(x){
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            crisprNuclease=SpCas9,
#     #                            bsgenome=bsgenome_human,
#     #                            n_mismatches=x)
#     # expect_true(max(mcols(out)$n_mismatches) <= x)
#     # out <- addSpacerAlignments(gs,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            n_mismatches=x)
#     # mm_cols <- paste0("n", 0:x)
#     # expect_true(all(mm_cols %in% colnames(mcols(out))))
#     # })
# })


# test_that("n_max_alignments arg must be a single positive integer value", {
#     ## argument passed on to crisprBowtie::runCrisprBowtie
#     # bad_input <- list(NA,
#     #                   -1,
#     #                   0,
#     #                   0.5)
#     #                   # c(1,1))  # allowed
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      crisprNuclease=SpCas9,
#     #                                      bsgenome=bsgenome_human,
#     #                                      n_max_alignments=x))
#     # })
# })


# test_that("all_alignments arg must be TRUE or FALSE", {
#     ## argument passed on to crisprBowtie::runCrisprBowtie
#     # bad_input <- list(NA,
#     #                   0)
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      crisprNuclease=SpCas9,
#     #                                      bsgenome=bsgenome_human,
#     #                                      all_alignments=x))
#     # })
# })


# test_that("crisprNuclease arg must be a CrisprNuclease object", {
#     # bad_input <- list(NA,
#     #                   "SpCas9")
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      crisprNuclease=x))
#     # })
# })


# test_that("crisprNuclease arg is appropriately applied", {
#     ## skipped due to long runtime
#     # out_spcas9 <- getSpacerAlignments(spacers=spacers,
#     #                                   aligner="bowtie",
#     #                                   aligner_index=bowtie_index,
#     #                                   bsgenome=bsgenome_human,
#     #                                   crisprNuclease=SpCas9)
#     # spcas9_motifs <- motifs(SpCas9, primary=TRUE, expand=TRUE,
#     #                         as.character=TRUE)
#     # expect_equal(unique(width(out_spcas9$spacer)), spacerLength(SpCas9))
#     # expect_equal(unique(width(out_spcas9$protospacer)), spacerLength(SpCas9))
#     # expect_true(all(as.character(out_spcas9$pam) %in% spcas9_motifs))
#     # 
#     # gr_input <- GRanges("chr12", IRanges(start=25224014, width=1000))
#     # gs <- findSpacers(gr_input, bsgenome=bsgenome_human, crisprNuclease=AsCas12a)
#     # out_ascas12a <- getSpacerAlignments(spacers=protospacers(gs),
#     #                                     aligner="bowtie",
#     #                                     aligner_index=bowtie_index,
#     #                                     bsgenome=bsgenome_human,
#     #                                     crisprNuclease=AsCas12a)
#     # ascas12a_motifs <- motifs(AsCas12a, primary=TRUE, expand=TRUE,
#     #                           as.character=TRUE)
#     # expect_equal(unique(width(out_ascas12a$spacer)), spacerLength(AsCas12a))
#     # expect_equal(unique(width(out_ascas12a$protospacer)), spacerLength(AsCas12a))
#     # expect_true(all(as.character(out_ascas12a$pam) %in% ascas12a_motifs))
# })


# test_that("canonical arg must be logical", {
#     # bad_input <- list(NULL,
#     #                   0,
#     #                   "TRUE")
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      crisprNuclease=SpCas9,
#     #                                      canonical=x))
#     # })
# })


# test_that("canonical arg is appropriately applied", {
#     ## skipped due to long runtime
#     # canonicalPams <- motifs(SpCas9, primary=TRUE, expand=TRUE,
#     #                         as.character=TRUE)
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            crisprNuclease=SpCas9,
#     #                            canonical=TRUE,
#     #                            n_mismatches=3)
#     # expect_true(all(as.character(out$pam) %in% canonicalPams))
#     # 
#     # noncanonicalPams <- motifs(SpCas9, primary=FALSE, expand=TRUE,
#     #                         as.character=TRUE)
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            crisprNuclease=SpCas9,
#     #                            canonical=FALSE,
#     #                            n_mismatches=3)
#     # expect_false(all(as.character(out$pam) %in% canonicalPams))
#     # expect_true(all(as.character(out$pam) %in% noncanonicalPams))
#     # 
#     # # ignore pams
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            crisprNuclease=SpCas9,
#     #                            canonical=NA,
#     #                            n_mismatches=3)
#     # expect_false(all(as.character(out$pam) %in% canonicalPams))
#     # expect_false(all(as.character(out$pam) %in% noncanonicalPams))
# })


# test_that("standard_chr_only arg must be TRUE or FALSE", {
#     ## skipped due to long runtime
#     # bad_input <- list(NULL,
#     #                   0,
#     #                   NA,
#     #                   "TRUE")
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      crisprNuclease=SpCas9,
#     #                                      standard_chr_only=x))
#     # })
# })


# test_that("standard_chr_only arg is appropriately applied", {
#     # standardChrs <- GenomeInfoDb::standardChromosomes(bsgenome_human)
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            crisprNuclease=SpCas9,
#     #                            standard_chr_only=FALSE)
#     # expect_false(all(seqlevels(out) %in% standardChrs))
# })


# test_that("both_strands arg must be TRUE or FALSE", {
#     # bad_input <- list(NULL,
#     #                   0,
#     #                   NA,
#     #                   "TRUE")
#     # lapply(bad_input, function(x){
#     #     expect_error(getSpacerAlignments(spacers=spacers,
#     #                                      aligner="biostrings",
#     #                                      custom_seq=custom_seq1,
#     #                                      crisprNuclease=SpCas9,
#     #                                      both_strands=x))
#     # })
# })


# test_that("both_strands arg is appropriately applied", {
#     # out_single_strand <- getSpacerAlignments(spacers=spacers,
#     #                                          aligner="biostrings",
#     #                                          custom_seq=custom_seq4,
#     #                                          crisprNuclease=SpCas9,
#     #                                          both_strands=FALSE)
#     # expect_true(length(out_single_strand) == 0)
#     # out_both_strands <- getSpacerAlignments(spacers=spacers,
#     #                                         aligner="biostrings",
#     #                                         custom_seq=custom_seq4,
#     #                                         crisprNuclease=SpCas9,
#     #                                         both_strands=TRUE)
#     # expect_true(length(out_both_strands) == length(spacers))
    
# })


# test_that("getSpacerAlignments returns a GRanges object", {
#     # out <- getSpacerAlignments(spacers=spacers,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            crisprNuclease=SpCas9)
#     # expect_true(is(out, "GRanges"))
# })


# test_that("getSpacerAlignments returns DNAStringSet columns for sequences", {
#     out <- getSpacerAlignments(spacers=spacers,
#                                aligner="bowtie",
#                                aligner_index=bowtie_index,
#                                bsgenome=bsgenome_human,
#                                crisprNuclease=SpCas9)
#     expect_true(is(mcols(out)$spacer, "DNAStringSet"))
#     expect_true(is(mcols(out)$protospacer, "DNAStringSet"))
#     expect_true(is(mcols(out)$pam, "DNAStringSet"))
# })


# test_that("getSpacerAlignments returns seqinfo compatible with bsgenome arg", {
#     out <- getSpacerAlignments(spacers=spacers,
#                                aligner="bowtie",
#                                aligner_index=bowtie_index,
#                                bsgenome=bsgenome_human,
#                                crisprNuclease=SpCas9)
#     expect_error(seqinfo <- checkCompatibleSeqinfo(out, bsgenome_human),
#                  regexp=NA)
#     expect_true(is(seqinfo, "Seqinfo"))
# })


# test_that("getSpacerAlignments returns non-null metadata", {
#     out <- getSpacerAlignments(spacers=spacers,
#                                aligner="bowtie",
#                                aligner_index=bowtie_index,
#                                bsgenome=bsgenome_human,
#                                crisprNuclease=SpCas9)
#     expect_true(!is.null(metadata(out)))
# })


# ## Specific to addSpacerAlignments and addSpacerAlignmentsIterative


# test_that("guideSet argument is required to be a GuideSet object", {
#     bad_input <- list("guideSetExample",
#                       as.data.frame(guideSetExample),
#                       as(guideSetExample, "GRanges"))
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(guideSet=x,
#                                          aligner="bowtie",
#                                          aligner_index=bowtie_index,
#                                          bsgenome=bsgenome_human))
#     })
# })


# test_that("an empty guideSet is handled gracefully", {
#     ## not yet implemented -- depends on calls to runCrisprBowtie, runCrisprBwa
#     # addSpacerAlignments(gs[0],
#     #                     aligner="bowtie",
#     #                     aligner_index=bowtie_index,
#     #                     bsgenome=bsgenome_human)
    
#     ## not yet implemented
#     # addSpacerAlignments(gs[0],
#     #                     aligner="biostrings",
#     #                     custom_seq=custom_seq1)
    
# })


# test_that("colname must be a character string", {
#     bad_input <- list(NULL,
#                       NA,
#                       0,
#                       TRUE,
#                       list(""),
#                       c("", ""))
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(gs,
#                                          aligner="biostrings",
#                                          custom_seq=custom_seq1,
#                                          colname=x))
#     })
# })


# test_that("colname is appended to mcols(guideSet) as GRangesList", {
#     colname <- "test_colname"
#     out <- addSpacerAlignments(gs,
#                                aligner="biostrings",
#                                custom_seq=custom_seq1,
#                                colname=colname)
#     expect_true(colname %in% colnames(mcols(out)))
#     expect_true(is(mcols(out)[[colname]], "GRangesList"))
# })


# test_that("addSummary controls addtion/updating of alignment summary columns", {
#     ## skipped due to long runtime
#     # sumCols <- paste0("n", 0:3)
#     # expect_false(any(sumCols %in% colnames(mcols(gs))))
#     # 
#     # out <- addSpacerAlignments(gs,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            addSummary=TRUE,
#     #                            n_mismatches=3)
#     # expect_true(all(sumCols %in% colnames(mcols(out))))
#     # alnSummary <- mcols(out)[, sumCols]
#     # out <- addSpacerAlignments(out,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index_mouse,
#     #                            bsgenome=bsgenome_mouse,
#     #                            addSummary=FALSE,
#     #                            n_mismatches=3)
#     # expect_equal(alnSummary, mcols(out)[, sumCols])
# })


# test_that("txObject arg must be a GRangesList or TxDb object", {
#     bad_input <- list(NA,
#                       grListExample[["transcripts"]])
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(gs,
#                                          aligner="bowtie",
#                                          aligner_index=bowtie_index,
#                                          bsgenome=bsgenome_human,
#                                          txObject=x))
#     })
# })


# test_that("nX_c columns added when txObject is provided", {
#     ## skipped due to long runtime
#     # out <- addSpacerAlignments(gs,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            txObject=grListExample,
#     #                            n_mismatches=3)
#     # cds_cols <- paste0("n", 0:3, "_c")
#     # expect_true(all(cds_cols %in% colnames(mcols(out))))
# })


# test_that("tssObject arg must be a valid GRanges object", {
#     bad_input <- list(NA,
#                       grListExample[["transcripts"]])
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(gs,
#                                          aligner="bowtie",
#                                          aligner_index=bowtie_index,
#                                          bsgenome=bsgenome_human,
#                                          tssObject=x))
#     })
# })


# test_that("nX_p columns added when tssObject is provided", {
#     ## skipped due to long runtime
#     # out <- addSpacerAlignments(gs,
#     #                            aligner="bowtie",
#     #                            aligner_index=bowtie_index,
#     #                            bsgenome=bsgenome_human,
#     #                            tssObject=tssObjectExample,
#     #                            n_mismatches=3)
#     # promoter_cols <- paste0("n", 0:3, "_p")
#     # expect_true(all(promoter_cols %in% colnames(mcols(out))))
# })


# test_that("anchor arg must be 'cut_site' or 'pam_site'", {
#     bad_input <- list(NA,
#                       0,
#                       character(0),
#                       # "cut",  # works -- partial matching
#                       "Cut_site")
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(gs,
#                                          aligner="bowtie",
#                                          aligner_index=bowtie_index,
#                                          bsgenome=bsgenome_human,
#                                          txObject=grListExample,
#                                          tssObject=tssObjectExample,
#                                          anchor=x))
#     })
# })


# test_that("anchor arg is appropriately applied", {
#     ## skipped due to long runtime
#     # out_pam_site <- addSpacerAlignments(gs,
#     #                                     aligner="bowtie",
#     #                                     aligner_index=bowtie_index,
#     #                                     bsgenome=bsgenome_human,
#     #                                     txObject=txdb_human,
#     #                                     anchor="pam_site",
#     #                                     n_mismatches=2)
#     # out_cut_site <- addSpacerAlignments(gs,
#     #                                     aligner="bowtie",
#     #                                     aligner_index=bowtie_index,
#     #                                     bsgenome=bsgenome_human,
#     #                                     txObject=txdb_human,
#     #                                     anchor="cut_site",
#     #                                     n_mismatches=2)
#     # out_pam_site <- alignments(out_pam_site)
#     # out_cut_site <- alignments(out_cut_site)
#     # expect_false(identical(out_pam_site$intergenic_distance,
#     #                        out_cut_site$intergenic_distance))
# })


# test_that("annotationType arg must be 'gene_symbol' or 'gene_id'", {
#     bad_input <- list(NA,
#                       0,
#                       character(0),
#                       "GENE_ID")
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignments(gs,
#                                          aligner="bowtie",
#                                          aligner_index=bowtie_index,
#                                          bsgenome=bsgenome_human,
#                                          txObject=grListExample,
#                                          tssObject=tssObjectExample,
#                                          annotationType=x))
#     })
# })


# test_that("annotationType arg is appropriately applied", {
#     # check alignmemnt colname output
# })


# test_that("tss_window arg must be an integer vector of length 2", {
#     ## skipped due to long runtime
#     # bad_input <- list(NA,
#     #                   0,
#     #                   c("1", "10"),
#     #                   list(1, 10),
#     #                   c(0.5, 10.5))
#     # lapply(bad_input, function(x){
#     #     expect_error(addSpacerAlignments(gs,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      tssObject=tssObjectExample,
#     #                                      tss_window=x))
#     # })
# })


# test_that("tss_window arg is appropriately applied", {
#     ## skipped due to long runtime
#     # out_no_hits <- addSpacerAlignments(gs,
#     #                                    aligner="bowtie",
#     #                                    aligner_index=bowtie_index,
#     #                                    bsgenome=bsgenome_human,
#     #                                    tssObject=tssObjectExample,
#     #                                    tss_window=c(0,0))
#     # expect_true(all(out_no_hits$n0_p == 0))
#     # expect_true(all(is.na(alignments(out_no_hits)$promoters)))
#     # 
#     # out_with_hits <- addSpacerAlignments(gs,
#     #                                      aligner="bowtie",
#     #                                      aligner_index=bowtie_index,
#     #                                      bsgenome=bsgenome_human,
#     #                                      tssObject=tssObjectExample,
#     #                                      tss_window=c(-500,500))
#     # expect_true(all(out_with_hits$n0_p > 0))
#     # expect_true(any(!is.na(alignments(out_with_hits)$promoters)))
# })


# test_that("alignmentThresholds must contain non-negative integers", {
#     bad_input <- list(NA,
#                       -1,
#                       "1",
#                       0.5,
#                       c(1,2),
#                       list(1))
#     lapply(bad_input, function(x){
#         expect_error(addSpacerAlignmentsIterative(gs,
#                                                   aligner="bowtie",
#                                                   aligner_index=bowtie_index,
#                                                   bsgenome=bsgenome_human,
#                                                   alignmentThresholds=c(n0=x)))
#     })
# })


# test_that("alignmentThresholds must not have duplicate names", {
#     # expect_error(addSpacerAlignmentsIterative(gs,
#     #                                           aligner="bowtie",
#     #                                           aligner_index=bowtie_index,
#     #                                           bsgenome=bsgenome_human,
#     #                                           alignmentThresholds=c(n0=1,
#     #                                                                 n1=2,
#     #                                                                 n1=3)))
# })


# test_that("alignmentThresholds are appropriately applied", {
#     ## skipped due to long runtime
#     # out <- addSpacerAlignmentsIterative(gs,
#     #                                     aligner="bowtie",
#     #                                     aligner_index=bowtie_index,
#     #                                     bsgenome=bsgenome_human,
#     #                                     n_mismatches=3,
#     #                                     alignmentThresholds=c(n0=0))
#     # expect_true(all(is.na(mcols(out)[["n1"]])))
#     # out <- addSpacerAlignmentsIterative(gs,
#     #                                     aligner="bowtie",
#     #                                     aligner_index=bowtie_index,
#     #                                     bsgenome=bsgenome_human,
#     #                                     n_mismatches=3,
#     #                                     alignmentThresholds=c(n0=1))
#     # good <- out$n0 <= 1
#     # expect_true(all(!is.na(mcols(out)$n1[good])))
#     # expect_true(all(is.na(mcols(out)$n1[!good])))
# })
