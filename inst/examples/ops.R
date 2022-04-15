library(crisprDesign)
library(crisprAnnotation)
lib <- getCrisprLibraryAnnotation("crisprko.cas9.mouse.gne5.wholegenome")
lib <- lib[,c("ID", "spacer_20mer", "rank", "gene_symbol")]
lib <- dplyr::rename(lib, spacer=spacer_20mer)
lib$opsBarcode <- substr(lib$spacer, 1,9)
ntcs <- lib[lib$gene_symbol=="NTC",]
ntcs$rank <- 1:nrow(ntcs)
genes <- unique(lib$gene_symbol)[1:100]
lib <- lib[lib$gene_symbol %in% genes,]
#lib <- lib[lib$rank<=4,]

opsLibrary <- designOpsLibrary(lib, gene_field="gene_symbol",)

# Adding ntcs:
opsLibrary <- updateOpsLibrary(opsLibrary,
                               df=ntcs,
                               n_guides=100,
                               gene_field="gene_symbol")
