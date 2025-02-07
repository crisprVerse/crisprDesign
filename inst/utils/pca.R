
#' @title Get principal components.
#' @description Get PCs from a Summarized Experiment object.
#' @param se A SummarizedExperiment object.
#' @param condition.field String specifying the name of the column in colData(se)
#'       that specifies the experimental groups.
#' @param reference.level String specifying the name of
#'     the reference group.
#' @param replicate.field String specifying the name of the replicate
#'     column in \code{colData(se)}.
#' @param how.many Numeric value indicating how many PCs to return.
#' 
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @export
getPCs <- function(se,
                   condition.field,
                   reference.level,
                   replicate.field,
                   how.many=5000
){
    how.many <- min(how.many, nrow(se))
    seRatios <- getLogRatios(se,
                             condition.field=condition.field,
                             reference.level=reference.level,
                             replicate.field=replicate.field)
    Y <- assays(seRatios)[[1]]
    # Scaling rows:
    Y   <- sweep(Y,1,rowMeans(Y, na.rm=TRUE),"-")
    wh  <- order(-rowVars(Y))[1:how.many]
    pca <- prcomp(Y[wh,,drop=FALSE],
                  scale=FALSE,
                  center=FALSE)
    pc1 <- pca$rotation[,1]
    pc2 <- pca$rotation[,2]
    percs <- round((pca$sdev^2)/sum(pca$sdev^2)*100,1)
    out <- list(pc1=pca$rotation[,1],
                pc2=pca$rotation[,2],
                perc1=percs[[1]],
                perc2=percs[[2]])
    return(out)
}
