library(crisprDesign)
library(crisprDesignGne)
library(crisprDesignData)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
txObject  <- txdb_human
tssObject <- tss_human
bsgenome  <- BSgenome.Hsapiens.UCSC.hg38
crisprNuclease <- crisprBase::SpCas9
bowtieIndex <- getBowtieIndex()
vcf <- system.file("extdata",
                   file="common_snps_dbsnp151_example.vcf.gz",
                   package="crisprDesign")


gr <- queryTxObject(txObject,
                    queryColumn="gene_symbol",
                    queryValue="KRAS",
                    featureType="cds")
start(gr) <- start(gr)-100


guides <- findSpacers(gr,
                      bsgenome=bsgenome,
                      crisprNuclease=crisprNuclease)
guides <- addSequenceFeatures(guides)
guides <- addRestrictionEnzymes(guides)
guides <- addSpacerAlignments(guides,
                              bsgenome=bsgenome,
                              aligner_index=bowtieIndex,
                              txObject=txObject,
                              tssObject=tssObject)
guides <- addOffTargetScores(guides)
guides <- addPamScores(guides)
guides <- addRepeats(guides, grRepeatsExample)
guides <- addSNPAnnotation(guides, vcf=vcf)
#guides <- addOnTargetScores(guides)
guides <- addLibraryUsage(guides)
guides <- addGeneAnnotation(guides,
                            txObject=txObject)
guides <- addTssAnnotation(guides,
                           tssObject=tssObject)









