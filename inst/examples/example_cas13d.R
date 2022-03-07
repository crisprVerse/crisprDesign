library(crisprDesignData)
library(crisprDesign)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
data(CasRx)
txObject <- txdb_human
txids <- c("ENST00000457313")
txids <- c("ENST00000367064") #CD55




out <- getMrnaSequences(txids,
                        bsgenome=bsgenome,
                        txObject=txObject)
guides <- findSpacers(out[1], crisprNuclease=CasRx)
#guides <- guides[1:10]
guides <- addSequenceFeatures(guides)
guides <- addPamScores(guides)
guides <- addRestrictionEnzymes(guides)


### Bowtie alignment:
bowtie_index="/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
guides <- addSpacerAlignments(guides,
                              txObject=txObject,
                              n_mismatches=1,
                              aligner="bowtie",
                              both_strand=TRUE,
                              aligner_index=bowtie_index)




aln <- alignments(guides)
n_mismatches=2
spacers <- spacers(guides, as.character=TRUE)


aln <- alignments(guides)








### Biostrings alignment:
load("../../../crisprDesignData/data/mrnasHuman.rda")
guides <- addSpacerAlignments(guides,
                              n_mismatches=2,
                              aligner="biostrings",
                              both_strand=TRUE,
                              custom_seq=mrnasHuman[1:10])








# On-target scoring:
guides <- addOnTargetScores(guides)


x <- 1:length(guides)
y <- guides$score_casrxrf
xlim=c(200,1400)
plot(x,y, pch=20, cex=0.5, col="grey55", xlim=xlim)
plot(x,y, pch=20, cex=0.5, col="grey55")
lines(predict(loess(y~x, span=0.03)), col=2, lwd=3)




