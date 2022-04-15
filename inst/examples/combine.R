library(crisprDesign)
library(crisprDesignData)
library(crisprDesignGne)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
data(tss_human)
txObject  <- txdb_human
tssObject <- tss_human
bsgenome  <- BSgenome.Hsapiens.UCSC.hg38
crisprNuclease <- crisprBase::SpCas9
index=getBowtieIndex()
vcf <- getSNPFile()
#vcf <- system.file("extdata",
#                   file="common_snps_dbsnp151_example.vcf.gz",
#                   package="crisprDesign")

gr1 <- queryTxObject(txObject,
                     queryColumn="gene_symbol",
                     queryValue="NRAS",
                     featureType="cds")
gr2 <- queryTxObject(txObject,
                     queryColumn="gene_symbol",
                     queryValue="KRAS",
                     featureType="cds")
grs <- list(gr1, gr2)


.getGuides <- function(gr){
    guides <- findSpacers(gr,
                          bsgenome=bsgenome,
                          crisprNuclease=crisprNuclease)
    guides <- addSequenceFeatures(guides)
    guides <- addRestrictionEnzymes(guides)
    guides <- addSpacerAlignments(guides,
                                  aligner_index=index,
                                  bsgenome=bsgenome,
                                  n_mismatches=1,
                                  txObject=txObject,
                                  tssObject=tssObject)
    guides <- addOffTargetScores(guides)
    guides <- addPamScores(guides)
    guides <- addRepeats(guides, grRepeatsExample)
    guides <- addSNPAnnotation(guides, vcf=vcf)
    guides <- addGeneAnnotation(guides,
                                txObject=txObject)
    guides <- addTssAnnotation(guides,
                               tssObject=tssObject)
    guides
}
grs <- lapply(grs, .getGuides)
old <- grs


combined <- Reduce("c", grs)












