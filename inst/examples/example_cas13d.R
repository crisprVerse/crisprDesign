library(crisprDesignData)
library(crisprDesign)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
data(CasRx)
txObject <- txdb_human
txids <- c("ENST00000457313")


out <- getMrnaSequences(txids,
                        bsgenome=bsgenome,
                        txObject=txObject)
guides <- findSpacers(out[1], crisprNuclease=CasRx)
#guides <- guides[1:10]
guides <- addOnTargetScores(guides)



load("../../../crisprDesignData/data/mrnasHuman.rda")
guides <- addSpacerAlignments(guides,
                              n_mismatches=2,
                              aligner="biostrings",
                              both_strand=TRUE,
                              custom_seq=mrnasHuman[1:10])

guides <- addSequenceFeatures(guides)
guides <- addPamScores(guides)
guides <- addRestrictionEnzymes(guides)




