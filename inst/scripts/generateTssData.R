library(crisprDesign)
library(GenomicRanges)
load("../../data/grListExample.rda")
tss <- grListExample[["fiveUTRs"]]
tss <- tss[tss$exon_rank==1]
cols <- c("tx_id", "gene_id", "gene_symbol")
mcols(tss) <- mcols(tss)[,cols]
tss$promoter <- paste0("P", seq_along(tss))
tss$ID <- paste0(tss$gene_symbol, "_", tss$promoter)
# Making sure we only retain one coordinate:
if (as.character(strand(tss))[[1]]=="+"){
    end(tss) <- start(tss)
} else {
    start(tss) <- end(tss)
}


tssObjectExample <- tss
save(tssObjectExample,
     file="../../data/tssObjectExample.rda")
