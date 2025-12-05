library(GenomicRanges)
library(GenomeInfoDb)
inputfile <- "/Users/fortin946/crisprVerseInputs/snps/gr.snps.db155.rds"
snps <- readRDS(inputfile)

load("../../data/grListExample.rda")
min <- min(start(grListExample$cds))-100
max <- max(end(grListExample$cds))+100
chr <- as.character(seqnames(grListExample$cds)[1])
gr <- GRanges(chr, IRanges(start=min, end=max))
genome(gr) <- "hg38"

snps <- snps[seqnames(snps)==chr]
snps <- snps[start(snps)>=min]
snps <- snps[end(snps)<max]
snpObjectExample <- snps
save(snpObjectExample, file="../../data/snpObjectExample.rda")
