data("SpCas9", "enAsCas12a", package="crisprBase")
bsgenome_human <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
bsgenome_mouse <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10

seq1 <- "CGAAATTTGCGTACCCATGTGACTGTGACGTCGATTGACTGCGTATGGCGGCATTTTCGATTC"
seq2 <- "CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTAGCCGGGTGTTA"
seq3 <- "DCATKACCCATDKTCADTTAYTCAGVTATGAANCATATCTAACACTAATVCCTTAGCTACGCATGGC"
seq4 <- "TTTACTGG"
seq5 <- ""

gr1 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(start=c(12345678, 12345678), width=500),
               strand=c("+", "-"))
gr2 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(start=c(12345678, 12345678), width=1),
               strand=c("+", "-"))
dnastringset1 <- DNAStringSet(c(seq1, seq2, seq3))
dnastring1 <- DNAString(seq1)


test_that("x must be a GRanges, DNAStringSet, DNAString or character vector", {
    bad_input <- list(NULL,
                      NA,
                      letters,
                      list(seq1),
                      data.frame(seq=seq1))
    lapply(bad_input, function(x){
        expect_error(findSpacers(x))
    })
    good_input <- list(gr1,
                       gr2,
                       dnastringset1,
                       dnastring1,
                       c(seq1=seq1, seq2=seq2, seq3=seq3))
    lapply(good_input, function(x){
        expect_error(findSpacers(x, bsgenome=bsgenome_human),
                     regexp=NA)
    })
})


test_that("crisprNuclease must be a crisprNuclease object or NULL", {
    bad_input <- list(NA,
                      "SpCas9")
    lapply(bad_input, function(x){
        expect_error(findSpacers(seq1, crisprNuclease=x))
    })
    good_input <- list(NULL,
                       SpCas9,
                       enAsCas12a)
    lapply(good_input, function(x){
        expect_error(findSpacers(seq1, crisprNuclease=x),
                     regexp=NA)
    })
})


test_that("bsgenome must be a BSgenome object or NULL when x is GRanges", {
    bad_input <- list(NA,
                      NULL,
                      "BSgenome.Hsapiens.UCSC.hg38")
    lapply(bad_input, function(x){
        expect_error(findSpacers(gr1, bsgenome=x))
        expect_error(findSpacers(seq1, bsgenome=x),
                     regexp=NA)
    })
    good_input <- list(bsgenome_human,
                       bsgenome_mouse)
    lapply(good_input, function(x){
        expect_error(findSpacers(gr1, bsgenome=x),
                     regexp=NA)
    })
})


test_that("if provided, bsgenome must match genome annotation for x", {
    genome(gr1) <- "hg38"
    expect_error(findSpacers(gr1, bsgenome=bsgenome_mouse))
    expect_error(findSpacers(gr1, bsgenome=bsgenome_human),
                 regexp=NA)
    expect_error(findSpacers(gr1, bsgenome=NULL))
})


test_that("spacer_len must be a single, positive numeric value or NULL", {
    bad_input <- list(NA,
                      -1,
                      0,
                      "1",
                      list(1),
                      c(1,1))
    lapply(bad_input, function(x){
        expect_error(findSpacers(seq1, spacer_len=x))
    })
    good_input <- list(NULL,
                       1,
                       10)
                       # 100) # gives error...
    lapply(good_input, function(x){
        expect_error(findSpacers(seq1, spacer_len=x),
                     regexp=NA)
    })
})


bad_booleans <- list(NULL, NA, "TRUE", 1)
good_booleans <- list(TRUE, FALSE)

test_that("canonical must be a single boolean value", {
    lapply(bad_booleans, function(x){
        expect_error(findSpacers(seq1, canonical=x))
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human, canonical=x))
    })
    lapply(good_booleans, function(x){
        expect_error(findSpacers(seq1, canonical=x),
                     regexp=NA)
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human, canonical=x),
                     regexp=NA)
    })
})

test_that("both_strands must be a single boolean value", {
    lapply(bad_booleans, function(x){
        expect_error(findSpacers(seq1, both_strands=x))
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human, both_strands=x))
    })
    lapply(good_booleans, function(x){
        expect_error(findSpacers(seq1, both_strands=x),
                     regexp=NA)
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human, both_strands=x),
                     regexp=NA)
    })
})

test_that("strict_overlap must be a single boolean value", {
    lapply(bad_booleans, function(x){
        expect_error(findSpacers(seq1, strict_overlap=x))
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human,
                                 strict_overlap=x))
    })
    lapply(good_booleans, function(x){
        expect_error(findSpacers(seq1, strict_overlap=x),
                     regexp=NA)
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human,
                                 strict_overlap=x),
                     regexp=NA)
    })
})

test_that("remove_ambiguities must be a single boolean value", {
    lapply(bad_booleans, function(x){
        expect_error(findSpacers(seq1, remove_ambiguities=x))
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human,
                                 remove_ambiguities=x))
    })
    lapply(good_booleans, function(x){
        expect_error(findSpacers(seq1, remove_ambiguities=x),
                     regexp=NA)
        expect_error(findSpacers(gr1, bsgenome=bsgenome_human,
                                 remove_ambiguities=x),
                     regexp=NA)
    })
})


test_that("x can be specified data structures", {
    expect_true(is(findSpacers(seq1), "GuideSet"))
    expect_true(is(findSpacers(seq2), "GuideSet"))
    expect_true(is(findSpacers(seq3), "GuideSet"))
    expect_true(is(findSpacers(seq4), "GuideSet"))
    expect_true(is(findSpacers(seq5), "GuideSet"))
    expect_true(is(findSpacers(gr1, bsgenome=bsgenome_human), "GuideSet"))
    expect_true(is(findSpacers(gr2, bsgenome=bsgenome_human), "GuideSet"))
    expect_true(is(findSpacers(dnastringset1), "GuideSet"))
    expect_true(is(findSpacers(dnastring1), "GuideSet"))
})


test_that("only GRanges input requires bsgenome arg", {
    expect_error(findSpacers(seq1),
                 regexp=NA)
    expect_error(findSpacers(seq1, bsgenome=bsgenome_human),
                 regexp=NA)
    expect_error(findSpacers(gr1)) # expect error
    expect_error(findSpacers(gr1, bsgenome=bsgenome_human),
                 regexp=NA)
    expect_error(findSpacers(dnastringset1),
                 regexp=NA)
    expect_error(findSpacers(dnastringset1, bsgenome=bsgenome_human),
                 regexp=NA)
    expect_error(findSpacers(dnastring1),
                 regexp=NA)
    expect_error(findSpacers(dnastring1, bsgenome=bsgenome_human),
                 regexp=NA)
})


test_that("function uses crisprNuclease argument appropriately", {
    default <- findSpacers(seq1, crisprNuclease=NULL)
    expect_equal(crisprNuclease(default), SpCas9)
    expect_true(all(as.character(default$pam) %in%
                        pams(SpCas9, as.character=TRUE)))
    expect_true(all(width(default$protospacer) == spacerLength(SpCas9)))
    
    spcas9 <- findSpacers(seq1, crisprNuclease=SpCas9)
    expect_identical(crisprNuclease(spcas9), SpCas9)
    expect_true(all(as.character(spcas9$pam) %in%
                        pams(SpCas9, as.character=TRUE)))
    expect_true(all(width(spcas9$protospacer) == spacerLength(SpCas9)))
    
    enascas12a <- findSpacers(seq1, crisprNuclease=enAsCas12a)
    expect_identical(crisprNuclease(enascas12a), enAsCas12a)
    expect_true(all(as.character(enascas12a$pam) %in%
                        pams(enAsCas12a, as.character=TRUE)))
    expect_true(all(width(enascas12a$protospacer) == spacerLength(enAsCas12a)))
})


test_that("function uses bsgenome argument appropriately", {
    human <- findSpacers(gr1, bsgenome=bsgenome_human)
    expect_equal(metadata(human)$bsgenome, bsgenome_human)
    human <- crisprBase::getProtospacerRanges(human, nuclease=SpCas9)
    protospacers <- paste0(protospacers(human, as.character=TRUE),
                           pams(human, as.character=TRUE))
    names(protospacers) <- names(human)
    expect_equal(getSeq(bsgenome_human, human, as.character=TRUE),
                 protospacers)
    
    mouse <- findSpacers(gr1, bsgenome=bsgenome_mouse)
    expect_equal(metadata(mouse)$bsgenome, bsgenome_mouse)
    mouse <- crisprBase::getProtospacerRanges(mouse, nuclease=SpCas9)
    protospacers <- paste0(protospacers(mouse, as.character=TRUE),
                           pams(mouse, as.character=TRUE))
    names(protospacers) <- names(mouse)
    expect_equal(getSeq(bsgenome_mouse, mouse, as.character=TRUE),
                 protospacers)
})


test_that("function uses spacer_len argument appropriately", {
    lens <- list(NULL,
                 1,
                 10,
                 100)
    lapply(lens, function(x){
        gs <- findSpacers(gr1,
                          crisprNuclease=SpCas9,
                          bsgenome=bsgenome_human,
                          spacer_len=x)
        expectedLength <- ifelse(is.null(x), spacerLength(SpCas9), x)
        expect_true(all(width(spacers(gs)) == expectedLength))
    })
})


test_that("function uses canonical argument appropriately", {
    canon <- findSpacers(gr1,
                         crisprNuclease=SpCas9,
                         bsgenome=bsgenome_human,
                         canonical=TRUE)
    noncanon <- findSpacers(gr1,
                            crisprNuclease=SpCas9,
                            bsgenome=bsgenome_human,
                            canonical=FALSE)
    expect_true(all(protospacers(canon, as.character=TRUE) %in%
                        protospacers(noncanon, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(noncanon, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(SpCas9, primary=TRUE, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(SpCas9, primary=FALSE, as.character=TRUE)))
    expect_true(all(pams(noncanon, as.character=TRUE) %in%
                        pams(SpCas9, primary=FALSE, as.character=TRUE)))
    expect_false(all(pams(noncanon, as.character=TRUE) %in%
                         pams(SpCas9, primary=TRUE, as.character=TRUE)))
    
    canon <- findSpacers(dnastringset1,
                         crisprNuclease=SpCas9,
                         canonical=TRUE)
    noncanon <- findSpacers(dnastringset1,
                            crisprNuclease=SpCas9,
                            canonical=FALSE)
    expect_true(all(protospacers(canon, as.character=TRUE) %in%
                        protospacers(noncanon, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(noncanon, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(SpCas9, primary=TRUE, as.character=TRUE)))
    expect_true(all(pams(canon, as.character=TRUE) %in%
                        pams(SpCas9, primary=FALSE, as.character=TRUE)))
    expect_true(all(pams(noncanon, as.character=TRUE) %in%
                        pams(SpCas9, primary=FALSE, as.character=TRUE)))
    expect_false(all(pams(noncanon, as.character=TRUE) %in%
                         pams(SpCas9, primary=TRUE, as.character=TRUE)))
})


test_that("function uses both_strands argument appropriately", {
    bothStrandSeq <- findSpacers(seq1, both_strands=TRUE)
    singleStrandSeq <- findSpacers(seq1, both_strands=FALSE)
    bothStrandGr <- findSpacers(gr1[1],
                                bsgenome=bsgenome_human,
                                both_strands=TRUE)
    singleStrandGr <- findSpacers(gr1[1],
                                  bsgenome=bsgenome_human,
                                  both_strands=FALSE)
    .getStrandType <- function(gs){
        strand <- as.character(strand(gs))
        unique(strand)
    }
    expect_true(all(c("-", "+") %in% .getStrandType(bothStrandSeq)))
    expect_true(all(c("+") %in% .getStrandType(singleStrandSeq)))
    expect_true(all(protospacers(singleStrandSeq, as.character=TRUE) %in%
                        protospacers(bothStrandSeq, as.character=TRUE)))
    expect_true(all(c("-", "+") %in% .getStrandType(bothStrandGr)))
    expect_true(all(c("+") %in% .getStrandType(singleStrandGr)))
    expect_true(all(protospacers(singleStrandGr, as.character=TRUE) %in%
                        protospacers(bothStrandGr, as.character=TRUE)))
})


test_that("function uses strict_overlap argument appropriately", {
    strict <- findSpacers(gr2, bsgenome=bsgenome_human,
                          strict_overlap=TRUE)
    not_strict <- findSpacers(gr2, bsgenome=bsgenome_human,
                              strict_overlap=FALSE)
    
    expect_equal(unique(start(gr2)), unique(strict$cut_site))
    expect_false(all(unique(start(gr2)) == unique(not_strict$cut_site)))
    expect_true(any(unique(start(gr2)) == unique(not_strict$cut_site)))
    overlaps <- findOverlaps(gr2,
                             getProtospacerRanges(not_strict, nuclease=SpCas9),
                             ignore.strand=TRUE)
    expect_equal(seq_along(not_strict), subjectHits(overlaps))
})


test_that("function uses remove_ambiguities argument appropriately", {
    remove <- findSpacers(seq3, remove_ambiguities=TRUE)
    keep <- findSpacers(seq3, remove_ambiguities=FALSE)
    
    expect_true(length(remove) == 0)
    expect_true(length(keep) > 0)
    ambiguousBases <- setdiff(Biostrings::DNA_ALPHABET, Biostrings::DNA_BASES)
    ambiguousBases <- intersect(ambiguousBases, LETTERS)
    pattern <- paste0("[", paste0(ambiguousBases, collapse=""), "]")
    expect_true(all(grepl(pattern, protospacers(keep))))
    
    gr <- GRanges(seqnames="chr1",
                  ranges=IRanges(start=9980, width=100),
                  strand="+")
    removeGr <- findSpacers(gr, crisprNuclease=enAsCas12a,
                            bsgenome=bsgenome_human, canonical=FALSE,
                            remove_ambiguities=TRUE)
    keepGr <- findSpacers(gr, crisprNuclease=enAsCas12a,
                          bsgenome=bsgenome_human, canonical=FALSE,
                          remove_ambiguities=FALSE)
    expect_true(length(keepGr) >= length(removeGr))
    expect_false(all(grepl(pattern, protospacers(removeGr))))
    expect_true(any(grepl(pattern, protospacers(keepGr))))
})


test_that("strict_overlap argument only applies to GRanges input", {
    strictGr <- findSpacers(gr2, bsgenome=bsgenome_human,
                            strict_overlap=TRUE)
    notStrictGr <- findSpacers(gr2, bsgenome=bsgenome_human,
                               strict_overlap=FALSE)
    strictSeq <- findSpacers(c(seq1, seq2, seq3), strict_overlap=TRUE)
    notStrictSeq <- findSpacers(c(seq1, seq2, seq3), strict_overlap=FALSE)
    
    expect_false(isTRUE(all.equal(strictGr, notStrictGr)))
    expect_equal(strictSeq, notStrictSeq)
})


test_that("strand=='*' is treated as both strands", {
    gr <- GRanges(seqnames="chr1",
                  ranges=IRanges(start=12345678, width=1000),
                  strand="*")
    bothStrand <- findSpacers(gr, bsgenome=bsgenome_human,
                              both_strands=TRUE)
    singleStrand <- findSpacers(gr, bsgenome=bsgenome_human,
                                both_strands=FALSE)
    expect_equal(bothStrand, singleStrand)
    expect_true(all(c("+", "-") %in% as.character(strand(singleStrand))))
})


test_that("input having no spacers is handled appropriately", {
    out <- findSpacers(seq4)
    expect_true(length(out) == 0)
    expect_true(is(out, "GuideSet"))
    
    out <- findSpacers(seq5)
    expect_true(length(out) == 0)
    expect_true(is(out, "GuideSet"))
    
    gr <- GRanges(seqnames="chr1",
                  ranges=IRanges(start=100, width=1),
                  strand="+")
    out <- findSpacers(gr, bsgenome=bsgenome_human)
    expect_true(length(out) == 0)
    expect_true(is(out, "GuideSet"))
    
})


test_that("GRanges are not extended outside acceptable chromosome limits", {
    gr <- GRanges(seqnames="chr1",
                  ranges=IRanges(start=1, width=1),
                  strand="+")
    expect_warning(out <- findSpacers(gr, bsgenome=bsgenome_human))
    expect_true(is(out, "GuideSet"))
})


test_that("provided names for x must be unique", {
    names(gr1) <- c('a', 'a')
    expect_error(findSpacers(gr1, bsgenome=bsgenome_human))
})
