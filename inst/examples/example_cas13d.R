library(crisprDesignData)
library(crisprDesign)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
data(CasRx)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
txObject <- txdb_human
txids <- c("ENST00000367064") #CD55


out <- getMrnaSequences(txids,
                        bsgenome=bsgenome,
                        txObject=txObject)
guides <- findSpacers(out[1], crisprNuclease=CasRx)
a=flattenGuideSet(guides)





guides <- addOnTargetScores(guides)



guideSet <- guides





#guides <- guides[1:10]
guides <- addSequenceFeatures(guides)
guides <- addPamScores(guides)
guides <- addRestrictionEnzymes(guides)


### Bowtie alignment:
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
guides <- addSpacerAlignments(guides,
                              txObject=txObject,
                              n_mismatches=3,
                              aligner="bowtie",
                              both_strand=TRUE,
                              aligner_index=bowtie_index)



# On-target scoring:
guides <- addOnTargetScores(guides)



### Biostrings alignment:
load("../../../crisprIndices/transcriptomes/human/mrnasHuman.rda")
guides <- addSpacerAlignments(guides,
                              n_mismatches=2,
                              aligner="biostrings",
                              both_strand=TRUE,
                              custom_seq=mrnasHuman[txids])

