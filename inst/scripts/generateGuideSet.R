library(crisprDesignS4)
library(GenomicRanges)
library(devtools)

# Getting guides in CDS of gene "IQSEC3"
gr <- queryTxObject(grListExample,
                    featureType="cds",
                    queryColumn="gene_symbol",
                    queryValue="IQSEC3")
genome(gr) <- "hg38"
guideSet <- findSpacers(gr)
guideSetExample <- guideSet
use_data(guideSetExample,
         compress="xz",
         overwrite=TRUE)
