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


test_that("empty ntcs argument returns unchanged object", {
    empty_ntcs <- all_ntcs[0]
    expect_error(addNtcs(guideSetExample, empty_ntcs),
                 regexp=NA)
})


test_that("ntcs are added to object and update the Seqinfo of object", {
    out <- addNtcs(guideSetExample, all_ntcs)
    expect_true(all(names(all_ntcs) %in% seqlevels(out)))
    expect_true(all(names(all_ntcs) %in% as.character(seqnames(out))))
    expect_true(all(names(all_ntcs) %in% names(out)))
    expect_true(!any(names(all_ntcs) %in% seqlevels(guideSetExample)))
    expect_true(!any(names(all_ntcs) %in% as.character(seqnames(guideSetExample))))
    expect_true(!any(names(all_ntcs) %in% names(guideSetExample)))
    extra_ntcs <- all_ntcs
    names(extra_ntcs) <- paste0("new_ntcs", seq_along(all_ntcs))
    extra_out <- addNtcs(out, extra_ntcs)
    expect_true(all(names(extra_ntcs) %in% seqlevels(extra_out)))
    expect_true(all(names(extra_ntcs) %in% as.character(seqnames(extra_out))))
    expect_true(all(names(extra_ntcs) %in% names(extra_out)))
    expect_true(!any(names(extra_ntcs) %in% seqlevels(out)))
    expect_true(!any(names(extra_ntcs) %in% as.character(seqnames(out))))
    expect_true(!any(names(extra_ntcs) %in% names(out)))
})


test_that("added ntcs have NA or empty list annotations", {
    out <- addNtcs(guideSetExampleFullAnnotation, all_ntcs)
    mcolnames <- colnames(mcols(out))
    listcols <- c("alignments", "enzymeAnnotation", "geneAnnotation",
                  "tssAnnotation", "snps")
    basecols <- c("protospacer", "pam", "pam_site")
    atomicCols <- setdiff(mcolnames, c(listcols, basecols))
    out_ntcs <- out[intersect(names(out), names(all_ntcs))]
    lapply(listcols, function(i){
        lapply(seq_along(out_ntcs), function(ii){
            if (is(mcols(out_ntcs)[[i]], "GRangesList")){
                expect_equal(length(mcols(out_ntcs)[[i]][[ii]]), 0)
            } else {
                expect_equal(nrow(mcols(out_ntcs)[[i]][[ii]]), 0)
            }
        })
    })
    lapply(atomicCols, function(x){
        expect_true(all(is.na(mcols(out_ntcs)[[x]])))
    })
})


test_that("ntcs names must be unique and distinct from object ids, seqnames", {
    intraduplicate <- all_ntcs
    names(intraduplicate)[2] <- names(intraduplicate)[1]
    interduplicate_id <- all_ntcs
    names(interduplicate_id)[1] <- names(guideSetExample)[1]
    interduplicate_seqname <- all_ntcs
    names(interduplicate_seqname) <- seqnames(guideSetExample)[1]
    bad_input <- list(intraduplicate,
                      interduplicate_id,
                      interduplicate_seqname)
    lapply(bad_input, function(x){
        expect_error(addNtcs(guideSetExample, x))
    })
    # adding same ntcs
    out <- addNtcs(guideSetExample, all_ntcs)
    expect_error(addNtcs(out, all_ntcs))
})



test_that("crisprDesign computational functions handle ntcs in GuideSet", {
    out <- addNtcs(guideSetExample, all_ntcs)
    expect_error(convertToMinMaxGRanges(out),
                 regex=NA)
    expect_error(convertToProtospacerGRanges(out),
                 regexp=NA) # no error, but meaningless output for ntcs
})






out <- addNtcs(head(guideSetExample), all_ntcs[1])
out_full <- addNtcs(head(guideSetExampleFullAnnotation), all_ntcs[1])


## ACCESSOR FUNCTIONS =========================================================

test_that("atomic accessor functions handle ntcs in GuideSet gracefully", {
    expect_error(spacers(out_full), regexp=NA)
    expect_error(pams(out_full), regexp=NA)
    expect_error(pamSites(out_full), regexp=NA)
    expect_error(cutSites(out_full), regexp=NA)
    expect_error(protospacers(out_full), regexp=NA)
})

test_that("snps accessor handles ntcs in GuideSet gracefully", {
    expect_error(snps(out_full),
                 regexp=NA)
    expect_false(is.null(snps(out_full)))
})

test_that("alignments accessors handles ntcs in GuideSet gracefully", {
    expect_error(alignments(out_full),
                 regexp=NA)
    expect_error(onTargets(out_full),
                 regexp=NA)
    expect_error(offTargets(out_full),
                 regexp=NA)
    expect_false(is.null(alignments(out_full)))
    expect_false(is.null(onTargets(out_full)))
    expect_false(is.null(offTargets(out_full)))
})

test_that("geneAnnotation accessor handles ntcs in GuideSet gracefully", {
    expect_error(geneAnnotation(out_full),
                 regexp=NA)
    expect_false(is.null(geneAnnotation(out_full)))
})

test_that("tssAnnotation accessor handles ntcs in GuideSet gracefully", {
    expect_error(tssAnnotation(out_full),
                 regexp=NA)
    expect_false(is.null(tssAnnotation(out_full)))
})

test_that("enzymeAnnotation accessor handles ntcs in GuideSet gracefully", {
    expect_error(enzymeAnnotation(out_full),
                 regexp=NA)
    expect_false(is.null(enzymeAnnotation(out_full)))
})

## guideSetExampleFullAnnotation lacks editedAlleles annotation
# test_that("editedAlleles accessor handles ntcs in GuideSet gracefully", {
#     expect_error(editedAlleles(out_full),
#                  regexp=NA)
#     expect_false(is.null(editedAlleles(out_full)))
# })



## ANNOTATION FUNCTIONS =======================================================

## split for each crisprDesign function (so more helpful message if test fails)
test_that("addCutSites handles ntcs in GuideSet gracefully", {
    expect_error(addCutSites(out),
                 regexp=NA)
})




test_that("addSNPAnnotation handles ntcs in GuideSet gracefully", {
    VCF_PATH <- system.file("extdata",
                            file="common_snps_dbsnp151_example.vcf.gz",
                            package="crisprDesign")
    expect_error(res <- addSNPAnnotation(out, vcf=VCF_PATH),
                 regexp=NA)
    expect_error(snps(res),
                 regexp=NA)
    expect_false(is.null(snps(res)))
})


test_that("addGeneAnnotation handles ntcs in GuideSet gracefully", {
    expect_error(res <- addGeneAnnotation(out, txObject=grListExample),
                 regexp=NA)
    expect_error(geneAnnotation(res),
                 regexp=NA)
    expect_false(is.null(geneAnnotation(res)))
})


test_that("addTssAnnotation handles ntcs in GuideSet gracefully", {
    expect_error(res <- addTssAnnotation(out, tssObject=tssObjectExample),
                 regexp=NA)
    expect_error(tssAnnotation(res),
                 regexp=NA)
    expect_false(is.null(tssAnnotation(res)))
})


test_that("addOffTargetScores handles ntcs in GuideSet gracefully", {
    expect_error(addOffTargetScores(out_full),
                 regexp=NA)
})


test_that("addPamScores handles ntcs in GuideSet gracefully", {
    expect_error(addPamScores(out),
                 regexp=NA)
})


test_that("addRestrictionEnzymes handles ntcs in GuideSet gracefully", {
    expect_error(res <- addRestrictionEnzymes(out),
                 regexp=NA)
    expect_error(enzymeAnnotation(res),
                 regexp=NA)
    expect_false(is.null(enzymeAnnotation(res)))
})


test_that("addSpacerAlignments/Iterative handles ntcs in GuideSet gracefully", {
    ## create bowtie index
    fasta  <- system.file(package="crisprDesign", "fasta/chr12.fa")
    outdir <- tempdir()
    Rbowtie::bowtie_build(fasta,
                          outdir=outdir,
                          force=TRUE,
                          prefix="tempIndex")
    index <- file.path(outdir, "tempIndex")
    expect_error(res <- addSpacerAlignments(
        out,
        txObject=grListExample,
        tssObject=tssObjectExample,
        aligner_index=index,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38),
                 regexp=NA)
    expect_error(addSpacerAlignmentsIterative(
        out,
        txObject=grListExample,
        tssObject=tssObjectExample,
        aligner_index=index,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38),
                 regexp=NA)
    expect_error(alignments(res),
                 regexp=NA)
    expect_false(is.null(alignments(res)))
})




test_that("addSequenceFeatures handles ntcs in GuideSet gracefully", {
    expect_error(addSequenceFeatures(out),
                 regexp=NA)
})


test_that("GuideSet with ntcs can be converted to data.frame", {
    expect_error(as.data.frame(out),
                 regexp=NA)
})


test_that("addEditingSites handles ntcs in GuideSet gracefully", {
    gs <- out
    metadata(gs)$CrisprNuclease <- BE4max
    expect_error(addEditingSites(gs, "C2T"),
                 regexp=NA)
})


test_that("addExonTable handles ntcs in the GuideSet gracefully", {
    outga <- addGeneAnnotation(out, txObject=grListExample)
    expect_error(addExonTable(outga,
                              gene_id="ENSG00000120645",
                              txObject=grListExample),
                 regexp=NA)
})


## uses local file
# test_that("addConservationScores handles ntcs in the GuideSet gracefully", {
#     conservationFile <- getConservationFiles("human")
#     expect_error(addConservationScores(out,
#                                        conservationFile=conservationFile),
#                  regexp=NA)
# })


test_that("addDistanceToTss handles ntcs in the GuideSet gracefully", {
    tss_id <- "ENSG00000120645_P1"
    expect_error(addDistanceToTss(out_full, tss_id),
                 regexp=NA)
})


test_that("addTxTable handles ntcs in the GuideSet gracefully", {
    gene_id <- "ENSG00000120645"
    expect_error(addTxTable(out_full, gene_id, grListExample),
                 regexp=NA)
})


## uses local files
test_that("addCrispraiScores handles ntcs in the GuideSet gracefully", {
    # gr <- queryTss(tssObjectExample,
    #                "gene_symbol",
    #                "IQSEC3")
    # gs <- findSpacers(gr,
    #                   crisprNuclease=SpCas9,
    #                   bsgenome=BSgenome.Hsapiens.UCSC.hg38)
    # gs <- addNtcs(head(gs), all_ntcs[1])
    # chromatinFiles <- "~/crisprIndices/chromatin/hg38"
    # chromatinFiles <- file.path(chromatinFiles, list.files(chromatinFiles))
    # names(chromatinFiles) <- c("dnase", "faire", "mnase")
    # fastaFile <- "~/crisprIndices/genomes/hg38/hg38.fa.gz"
    # addCrispraiScores(gs,
    #                   gr=gr,
    #                   tssObject=tssObjectExample,
    #                   chromatinFiles=chromatinFiles,
    #                   fastaFile=fastaFile)
})


test_that("addEditedAlleles handles ntcs in the GuideSet gracefully", {
    gs <- out
    metadata(gs)$CrisprNuclease <- BE4max
    txTable <- getTxInfoDataFrame(tx_id="ENST00000538872",
                                  txObject=grListExample,
                                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)
    
    expect_error(res <- addEditedAlleles(gs,
                                         baseEditor=BE4max,
                                         txTable=txTable,
                                         editingWindow=c(-20, -8)),
                 regexp=NA)
    expect_error(editedAlleles(res),
                 regexp=NA)
    expect_false(is.null(editedAlleles(res)))
})


test_that("addIsoformAnnotation handles ntcs in the GuideSet gracefully", {
    tx_id <- "ENST00000538872"
    expect_error(addIsoformAnnotation(out_full,
                                      tx_id="ENST00000538872"),
                 regexp=NA)
})


test_that("addPfamDomains handles ntcs in the GuideSet gracefully", {
    # pfamTable <- preparePfamTable(grListExample,
    #                               mart_dataset="hsapiens_gene_ensembl")
    # expect_error(addPfamDomains(out_full,
    #                             pfamTable=pfamTable),
    #              regexp=NA)
})

