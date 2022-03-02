library(crisprDesignData)
library(crisprDesign)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

txObject <- txdb_human
txids <- c("ENST00000457313")
out <- getMrnaSequences(txids,
                        bsgenome=bsgenome,
                        txObject=txObject)
data(CasRx)
data(SpCas9)

guides <- findSpacers(out[1],
                      crisprNuclease=CasRx, spacer_len=23)
guides <- addSequenceFeatures(guides)
guides <- guides[1]


load("../../../crisprDesignData/data/mrnasHuman.rda")

guides <- addSpacerAlignments(guides,
                              n_mismatches=0,
                              aligner="biostrings",
                              both_strand=TRUE,
                              custom_seq=mrnasHuman)





