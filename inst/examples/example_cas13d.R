library(crisprDesignData)
library(crisprDesign)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

txObject <- txdb_human
txids <- c("ENST00000457313")
out <- crisprDesign:::getMrnaSequences(txids,
     bsgenome=bsgenome,
     txObject=txObject)
data(CasRx)
data(SpCas9)

guides <- findSpacers(out[1], crisprNuclease=CasRx, spacer_len=23)

#guides <- findSpacers(out[1], crisprNuclease=SpCas9, spacer_len=23)
spacers(guides)
protospacers(guides)
