library(crisprDesign)
library(GenomicRanges)
load("../../data/grListExample.rda")
tss <- grListExample[["fiveUTRs"]]
tss <- tss[tss$exon_rank==1]
cols <- c("tx_id", "gene_id", "gene_symbol")
mcols(tss) <- mcols(tss)[,cols]
tss$ID <- paste0(tss$gene_symbol, "_P", seq_along(tss))
# Making sure we only retain one coordinate:
if (as.character(strand(tss))[[1]]=="+"){
    end(tss) <- start(tss)
} else {
    start(tss) <- end(tss)
}


tssObjectExample <- tss
save(tssObjectExample,
     file="../../data/tssObjectExample.rda")
