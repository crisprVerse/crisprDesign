#' @export
equalizeSkoConstructs <- function(se,
                                   type.field="class",
                                   type.levels.sko1="target_neg",
                                   type.levels.sko2="neg_target",
                                   gene1.field="gene_symbol_1",
                                   gene2.field="gene_symbol_2",
                                   gene.pair.field="group",
                                   sko.control="neg"
){

    ann <- rowData(se)
    ann1 <- ann[ann[[type.field]] %in% type.levels.sko1,,drop=FALSE]
    ann2 <- ann[ann[[type.field]] %in% type.levels.sko2,,drop=FALSE]
    commonGuides <- intersect(ann1[["ID_1"]], ann2[["ID_2"]])
    if (length(commonGuides)==0){
        stop("No common guides at position 1 and position 2.")
    }
    id1 <- ann1$ID[match(commonGuides, ann1[["ID_1"]])]
    id2 <- ann2$ID[match(commonGuides, ann2[["ID_2"]])]
    wh1 <- match(id1, ann$ID)
    wh2 <- match(id2, ann$ID)

    # Normalization:
    assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    Y <- assays(se)[[1]]
    Y <- log2(Y+1)
    for (col in seq_len(ncol(Y))){
        target <- (Y[wh1,col]+Y[wh2,col])/2
        Y[wh1,col] <- target
        Y[wh2,col] <- target
    }
    Y <- 2^Y-1
    assays(se)[[1]] <- Y
    return(se)
}



# type.field="class"
# type.levels.dko="target_target"
# gene1.field="gene_symbol_1"
# gene2.field="gene_symbol_2"
# gene.pair.field="group"
# sko.control="neg"
# min.pairs=60
# min.per.position=5

#' @export
equalizeDkoConstructs <- function(se,
                                  type.field="class",
                                  type.levels.dko="target_target",
                                  gene1.field="gene_symbol_1",
                                  gene2.field="gene_symbol_2",
                                  gene.pair.field="group",
                                  sko.control="neg",
                                  min.pairs=60,
                                  min.per.position=5
){

    ann <- rowData(se)
    ann.dko <- ann[ann[[type.field]] %in% type.levels.dko,]
    ann.dko <- ann.dko[, c(gene1.field, gene2.field)]
    ann.dko <- ann.dko[!duplicated(ann.dko),]
    genes <- union(ann.dko[[gene1.field]], ann.dko[[gene2.field]])
    inPos1 <- genes %in% ann.dko[[gene1.field]]
    inPos2 <- genes %in% ann.dko[[gene2.field]]
    genesInBoth <- genes[inPos1 & inPos2]
    ns <- vapply(genesInBoth, function(x){
        sum(ann.dko[[gene1.field]]==x | ann.dko[[gene2.field]]==x)
    }, FUN.VALUE=0)
    genes <- genesInBoth[ns>=min.pairs]
    ns1 <- vapply(genes, function(x){
        sum(ann.dko[[gene1.field]]==x)
    }, FUN.VALUE=0)
    ns2 <- vapply(genes, function(x){
        sum(ann.dko[[gene2.field]]==x)
    }, FUN.VALUE=0)
    good <- ns1>=min.per.position & ns2>=min.per.position
    genes <- sort(genes[good])

    # Subsetting object we need to work on:
    ann <- rowData(se)
    cond1 <- ann[[type.field]] %in% type.levels.dko
    cond2 <- ann[[gene1.field]] %in% genes
    cond3 <- ann[[gene2.field]] %in% genes
    wh.subset <- cond1 & (cond2|cond3)
    see <- se[wh.subset,]
    ann <- rowData(see)

    # Normalization:
    Y <- as.matrix(assays(see)[[1]])
    Y <- log2(Y+1)
    for (gene in genes){
        wh1 <- which(ann[[gene1.field]]==gene)
        wh2 <- which(ann[[gene2.field]]==gene)
        for (col in seq_len(ncol(Y))){
            c1 <- median(Y[wh1,col])
            c2 <- median(Y[wh2,col])
            cc <- mean(c(c1,c2))
            c1 <- c1 - cc
            c2 <- c2 - cc
            Y[wh1,col] <- Y[wh1,col] - c1
            Y[wh2,col] <- Y[wh2,col] - c2
        }
    }
    Y <- 2^Y-1
    assays(se)[[1]][wh.subset,] <- Y
    return(se)
}


