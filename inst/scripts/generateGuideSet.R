library(crisprDesign)
library(GenomicRanges)
library(devtools)
data(grListExample)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38


# Getting guides in CDS of gene "IQSEC3"
gr <- queryTxObject(grListExample,
                    featureType="cds",
                    queryColumn="gene_symbol",
                    queryValue="IQSEC3")
guideSet <- findSpacers(gr,
                        bsgenome=bsgenome)

guideSetExample <- guideSet
use_data(guideSetExample,
         compress="xz",
         overwrite=TRUE)
