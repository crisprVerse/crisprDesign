library(crisprDesignS4)
library(GenomicRanges)
txdb <- getTxDb(organism="Homo sapiens")
grList <- TxDb2GRangesList(txdb)
meta <- metadata(grList)
#gene <- "ENSG00000133703"
gene="ENSG00000120645" #IQSEC3
grList <- lapply(grList, function(gr){
    gr[gr$gene_id==gene]
})
grListExample <- GRangesList(grList)
metadata(grListExample) <- meta
GenomeInfoDb::genome(grListExample) <- "hg38"
save(grListExample,
     file="../../data/grListExample.rda")
