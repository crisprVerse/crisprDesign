library(crisprDesignData)
library(GenomicRanges)
gr <- gr.repeats.hg38
load("../../data/grListExample.rda")
hits <- queryHits(findOverlaps(gr, grListExample))
hits <- unique(hits)
gr <- gr[hits]
grRepeatsExample <- gr
save(grRepeatsExample,
     file="../../data/grRepeatsExample.rda")
