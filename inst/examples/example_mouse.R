library(crisprDesignS4)
library(crisprDesignDataS4)
library(GenomeInfoDb)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10 

gr <- queryTxObject(txdb_mouse,
                    queryColumn="gene_symbol",
                    queryValue="Kras")
GenomeInfoDb::genome(gr) <- "mm10"
guides <- findSpacers(gr,
                      bsgenome=bsgenome)
seqlengths(gr)[1:10]
seqlengths(bsgenome)[1:10]