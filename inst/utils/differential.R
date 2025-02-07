#' @title Get single sample Rho statistic
#' @description Get single sample Rho statistic.
#' 
#' @param lfc A numeric vector of log-fold changes, or any other continuous statistics
#' @param annotation CRISPR library annotation. Number of rows must equal length of lfc
#' @param permutations Number of permutations to calculate p-values
#' @param permutation.seed Permutation seed for reproducibility
#' 
#' @return A data.frame containing medianlogFC, and Rho and P-values for the enrichment
#' and depletion tests. 
#' @importFrom gCrisprTools ct.prepareAnnotation
#' @importFrom gCrisprTools ct.applyAlpha
#' @importFrom gCrisprTools ct.RRAaPvals
#' @importFrom gCrisprTools ct.RRAalpha
#' @importFrom stats p.adjust
#' @export
getSSRhoStatistic <- function(lfc,
                              annotation = NULL,
                              permutations=1000,
                              permutation.seed = NULL
){
  annotation <- as.data.frame(annotation)
  cols <- colnames(annotation)
  if (!"geneSymbol" %in% cols){
    if ("gene_symbol"  %in% cols){
      annotation$geneSymbol=annotation[, "gene_symbol"]
    } else {
      annotation$geneSymbol=annotation[, "gene.symbol"]
    }
  }
  annotation$geneID <- annotation$geneSymbol
  annotation$geneID[annotation$geneSymbol=="NTC"] <- "no_gid"
  stats <- cbind(lfc, lfc, lfc)
  key <- gCrisprTools::ct.prepareAnnotation(annotation)
  key <- data.frame(key)
  rra.input <- gCrisprTools::ct.applyAlpha(stats, RRAalphaCutoff = 0.1, scoring ="fc")
  geneP.depletion <- gCrisprTools::ct.RRAaPvals(rra.input[, "scores.deplete", 
          drop = FALSE], 
          g.key = key, 
          permute = permutations, 
          permutation.seed = permutation.seed
  )
  geneP.enrichment <- gCrisprTools::ct.RRAaPvals(rra.input[, "scores.enrich", 
          drop = FALSE],
          g.key = key,
          permute = permutations,
          permutation.seed = permutation.seed
  )
  rhoEnrich <- gCrisprTools::ct.RRAalpha(rra.input[, "scores.enrich", drop = FALSE], 
      g.key = key,
      shuffle = FALSE
  )
  rhoDeplete <- gCrisprTools::ct.RRAalpha(rra.input[, "scores.deplete", drop = FALSE], 
      g.key = key,
      shuffle = FALSE
  )
  # Let's add logFC:
  lfcs <- split(lfc, annotation$geneSymbol)
  lfcs <- sapply(lfcs, median, na.rm=TRUE)
  lfcs <- lfcs[match(names(lfcs), names(rhoEnrich))]
  fdr.enrichment <- p.adjust(geneP.enrichment, method="fdr")
  fdr.depletion  <- p.adjust(geneP.depletion, method="fdr")
  names(lfcs)    <- names(rhoDeplete)
  data.frame(
    median.lfc=lfcs,
    rho.depletion=rhoDeplete,
    pval.depletion=geneP.depletion,
    rho.enrichment=rhoEnrich,
    pval.enrichment=geneP.enrichment,
    fdr.enrichment=fdr.enrichment,
    fdr.depletion=fdr.depletion)
}





#' @title Add a standardized timepoint column in SE 
#' @description Add a standardized timepoint column in SE
#' @param se A SummarizedExperiment object.
#' @return A SummarizedExperiment with a "Timepoint" column added to colData 
#' @import SummarizedExperiment
#' @importFrom magrittr %>%
#' @export
createTimepoint <- function(se){
    if("DelayedMatrix" %in% class(assays(se)[[1]])){
      assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    } 
    if (!"Day" %in% colnames(colData(se))){
      stop("Day should be provided in colData(se)")
    }
    day <- colData(se)$Day 
    day <- gsub("Day|day|D", "", day) %>% as.numeric
    if (!"Replicate" %in% colnames(colData(se))){
      stop("Replicate should be provided in colData(se)")
    }
    replicate  <- colData(se)$Replicate %>% as.character
    replicates <- unique(replicate)
    time <- day
    for (i in 1:length(replicates)) {
        indices <- which(replicate == replicates[[i]])
        temp <- factor(time[indices], ordered = TRUE)
        levels(temp) <- 1:length(levels(temp))
        time[indices] <- temp
    }
    colData(se)$Timepoint <- time
    se
}





