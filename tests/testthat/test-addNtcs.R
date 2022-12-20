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
            expect_equal(length(mcols(out_ntcs)[[i]][[ii]]), 0)
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
    # getPAMSequence(seqnames(out), pamSites(out), strand(out)) # error
    # getSpacerSequence(seqnames(out), pamSites(out), strand(out)) # error
    # convertToMinMaxGRanges(out) # issue with strand
    expect_error(convertToProtospacerGRanges(out),
                 regexp=NA) # no error, but meaningless output for ntcs
})





data("txdb_human", package="crisprDesignData")
data("tss_human", package="crisprDesignData")

out <- addNtcs(head(guideSetExample), all_ntcs[1])
out_full <- addNtcs(head(guideSetExampleFullAnnotation), all_ntcs[1])




## split for each crisprDesign function (so more helpful message if test fails)
test_that("addCutSites handles ntcs in GuideSet gracefully", {
    expect_error(addCutSites(out),
                 regexp=NA)
})


test_that("addRepeats handles ntcs in GuideSet gracefully", {
    data("gr.repeats.hg38", package="crisprDesignData")
    expect_error(addRepeats(out, gr.repeats=gr.repeats.hg38),
                 regexp=NA)
})


test_that("addSNPAnnotation handles ntcs in GuideSet gracefully", {
    VCF_PATH <- system.file("extdata",
                            file="common_snps_dbsnp151_example.vcf.gz",
                            package="crisprDesign")
    expect_error(addSNPAnnotation(out, vcf=VCF_PATH),
                 regexp=NA)
})


test_that("addGeneAnnotation handles ntcs in GuideSet gracefully", {
    expect_error(addGeneAnnotation(out, txObject=txdb_human),
                 regexp=NA)
})


test_that("addTssAnnotation handles ntcs in GuideSet gracefully", {
    expect_error(addTssAnnotation(out, tssObject=tss_human),
                 regexp=NA)
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
    expect_error(addRestrictionEnzymes(out),
                 regexp=NA)
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
    expect_error(addSpacerAlignments(
        out,
        txObject=txdb_human,
        tssObject=tss_human,
        aligner_index=index,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38),
                 regexp=NA)
    expect_error(addSpacerAlignmentsIterative(
        out,
        txObject=txdb_human,
        tssObject=tss_human,
        aligner_index=index,
        bsgenome=BSgenome.Hsapiens.UCSC.hg38),
                 regexp=NA)
})


test_that("addOnTargetScores handles ntcs in GuideSet gracefully", {
    expect_error(addOnTargetScores(out, methods=c("deephf")),
                 regexp=NA) # need to test all methods?
})


test_that("addSequenceFeatures handles ntcs in GuideSet gracefully", {
    expect_error(addSequenceFeatures(out),
                 regexp=NA)
})


test_that("GuideSet with ntcs can be converted to data.frame", {
    expect_error(as.data.frame(out),
                 regexp=NA)
})





## other crisprDesign functions that may need testing
# addEditingSites
# addExonTable_consensusIsoform
# addConservationScores
# addDistanceToTss
# addExonTable
# addTxTable
# addCrispraiScores
# addEditedAlleles
# addExonTable_allIsoforms
# addIsoformAnnotation
# addPfamDomains

## check that function can also handle PairedGuideSet?
