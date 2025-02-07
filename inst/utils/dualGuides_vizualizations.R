#' @title Generate Boxplot with scatterplot overlay of two
#'     given genes of interest from dualGuideDLFC logr output
#' @description Generate Boxplot/scatterplot of two given genes of
#'     interest from dualGuideDLFC logr output.
#' @param se Either a numeric matrix or a SummarizedExperiment.
#' @param gene.pair List of character vector of gene pair of interest
#'     (e.g. list(c("G1", "G2"))).
#' @param xlim should be number of boxplots (guides) to be plotted.
#' @param ylim y axis limits, the function gathers this automatically.
#' @param aspect.ratio default is 2/4, which is suitable for most plots.
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_jitter position_jitter 
#' @importFrom ggplot2 theme theme_bw element_text geom_vline geom_hline ggtitle labs ylim
#' @importFrom stats relevel
#' @export
triplot <- function(se,
                    gene1,
                    gene2,
                    gene.pair=NULL,
                    sko.control="neg",
                    dko.naive=NULL,
                    dko.bliss=NULL,
                    assay=1,
                    xlim=NULL,
                    ylim=NULL,
                    aspect.ratio=9/5,
                    main="",
                    use.shape=TRUE){
    if (!is.null(gene.pair)){
        x <- strsplit(gene.pair, split="_")[[1]]
        gene1 <- x[[1]]
        gene2 <- x[[2]]
    }
    genes1 <- listGenes(se, sko.position="first")
    genes2 <- listGenes(se, sko.position="second")
    if ((!gene1 %in% genes1) | (!gene2 %in% genes2)){
        temp  <- gene1
        gene1 <- gene2
        gene2 <- temp
    }
    if ((!gene1 %in% genes1) | (!gene2 %in% genes2)){
        stop("At least one gene does not exist in the annotation.")
    }
    sko1 <- .getSkoDataBX(se=se, gene=gene1, assay=assay, sko.position="first")
    sko2 <- .getSkoDataBX(se=se, gene=gene2, assay=assay, sko.position="second")
    dko  <- .getDkoDataBX(se=se, gene1=gene1, gene2=gene2, assay=assay)
    all.data <- rbind(sko1, dko, sko2)
    all.data$ID <- rownames(all.data)
    all.data$Guide <- gsub("_[0-9]+|_[a-z]+$", "", all.data$Label)
    all.data$Guide <- factor(all.data$Guide, levels=unique(all.data$Guide))

    if (use.shape){
        # Adding information about guide1 and guide2:
        ann <- rowData(se)
        ann <- ann[match(all.data$ID, ann$ID),]
        all.data$ID_1 <- ann$ID_1
        all.data$ID_2 <- ann$ID_2
        col <- as.factor(all.data$ID_1)
        col <- relevel(col, ref=sko.control)
        col <- as.numeric(col)

        # Giving a shape based on sgRNA2:
        shapes <- c(16,15,17,8,11,9)
        shapes <- c(shapes, setdiff(1:25, shapes))
        shape  <- as.factor(all.data$ID_2)
        shape  <- relevel(shape, ref=sko.control)
        shape  <- shapes[as.numeric(shape)]
    } else {
        col   <- as.numeric(as.factor(all.data$Guide)) # color based on guide
        shape <- rep(20, length(col))
    }

    

    p <- ggplot(all.data, aes_string(x="Guide",
                                     y="logRatios",
                                     Guide=FALSE)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(position=position_jitter(width=.1, height=0),
                    color=col,
                    shape=shape) +
        theme_bw() +
        theme(aspect.ratio=aspect.ratio,
              axis.text.x=element_text(angle=90,
                                       vjust=0.5,
                                       hjust=1)) +
        geom_vline(xintercept=1.5) +
        geom_vline(xintercept=2.5) +
        geom_hline(yintercept=0) + 
        ggtitle(main) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(x = "") 
    if (!is.null(xlim)){
        p <- p + xlim(xlim[1], xlim[2])
    }
    if (!is.null(ylim)){
        p <- p + ylim(ylim[1], ylim[2])
    }
    if (!is.null(dko.naive)){
        p <- p + geom_hline(yintercept=dko.naive, colour="orange") 
    }
    if (!is.null(dko.bliss)){
        p <- p + geom_hline(yintercept=dko.bliss, colour="red") 
    }
    return(p)
}





#' @title Generates boxplot of specified classes of guides between two samples.
#' @description Generates boxplot of specified classes of guides between two samples.
#' @param se.ratios A Summarized Experiment object or matrix.
#' @param class.field Character string indicating column in rowData(se.ratios) that specifies class.
#' @param ylab Character string specifying Y axis label.
#' @param classes Character string indicating which classes to plot.
#' @param ... Arguments to be passed to \code{boxplot}.
#' @return Returns nothing. a plot is produced as a side effect.
#' @export
#' @importFrom graphics boxplot
drawDualGuideBoxplotsQC <- function(se.ratios,
                                    class.field="class",
                                    ylab="Log2 fold change", 
                                    classes=unique(rowData(se.ratios)[[class.field]]),
                                    ...
){
    ann <- rowData(se.ratios)
    if (is(se.ratios, "SummarizedExperiment")) {
        if ("logRatios" %in% names(assays(se.ratios))) {
            Y <- assays(se.ratios)[["logRatios"]]
        }
        else if ("counts" %in% names(assays(se.ratios))) {
            Y <- assays(se.ratios)[["counts"]]
        }
        else {
            Y <- assays(se.ratios)[[1]]
        }
    }
    if (!all(classes %in% ann[[class.field]])){
        stop("classes not found in class.field column of rowData(se)")
    }
    xs <- split(Y, f = ann[[class.field]])[classes, drop = FALSE]
    boxplot(xs,las=2, outline=FALSE, ylab=ylab, ...)
    abline(h = 0)
}




#' @title Get class Key.
#' @description Gets colors for dual guide boxplot qc.
#' @param object List of classes with associated counts.
.getClassKey <- function(object){
    col     <- c("grey54", "grey0", "skyblue1", "skyblue2", "skyblue4", "lightpink", "lightpink1", "firebrick1", "darkolivegreen1", "darkolivegreen2", "darkolivegreen4")
    classes <- c("ntc_ntc", "neg_neg", "or_neg", "neg_or", "or_or", "neg_target", "target_neg", "target_target", "eg_neg", "neg_eg", "eg_eg")
    key     <- data.frame(colors=col,
                          classes=classes,
                          stringsAsFactors=FALSE)
    return(key)
}



#' @importFrom SummarizedExperiment colData
#' @export
drawDualQCBarplots <- function(se,
                               ...
){
    pheno <- as.data.frame(colData(se))
    col  <- 1:5
    if ("none" %in% colnames(pheno)){
        pheno$unmapped <- pheno$none
    }
    if ("nmapped" %in% colnames(pheno)){
        pheno$mapped <- pheno$nmapped
    }
    cols <- c("unmapped",
              "barcode1.only",
              "barcode2.only", 
              "invalid.pair", 
              "mapped")
    labels <- c("Unmapped", 
                "Valid barcode at position 1 only",
                "Valid barcode at position 2 only",
                "Invalid pair of barcodes",
                "Valid pair of barcodes")

    par(mar=c(10,7,8,15), xpd=TRUE)
    X <- as.matrix(pheno[,cols,drop=FALSE])
    X <- X/rowSums(X)*100
    ylim=c(0,100)
    barplot(t(X),
            ylab="",
            las=2,
            col=col,
            border=NA,
            ylim=ylim,
            ...)
    legend("topright",
           bty="n",
           legend=labels,
           col=col,
           fill=col,
           border=NA,
           inset=c(-0.6,0.5))
    mtext(side=2, "Number of reads", line=5)
}




.getSkoDataBX <- function(se,
                          gene,
                          sko.position=c("first", "second"),
                          assay=1,
                          ...){
    sko.position <- match.arg(sko.position)
    se <- getSkoData(se=se,
                     genes=gene,
                     sko.position=sko.position,
                     ...)
    ann <- rowData(se)
    # get everything after last underscore
    ann$col <- gsub(".*\\_*\\_", "", ann$ID)
    if (sko.position=="first"){
        Guide <- ann$ID_1
    } else if (sko.position=="second"){
        Guide <- ann$ID_2
    }
    Y <- assays(se)[[assay]]
    y <- rowMeans(Y, na.rm=TRUE)
    df   <- data.frame(Label=Guide,
                       Guide=Guide,
                       logRatios=y,
                       col=ann$col)
    return(df)
}






.getDkoDataBX <- function(se,
                          gene1,
                          gene2,
                          assay=1,
                          ...
){
    se <- getDkoData(se=se,
                     gene.pairs=list(c(gene1, gene2)),
                     return.matrix=FALSE,
                     ...)
    ann   <- rowData(se)
    ann$col   <- gsub(".*\\_*\\_", "", ann$ID)
    Y <- assays(se)[[assay]]
    y <- rowMeans(Y, na.rm=TRUE)
    df <- data.frame(Label=ann$ID,
                     Guide=ann$group,
                     logRatios=y,
                     col=ann$col)
    return(df)
}


.getExpectedDko <- function(ratios,
                            assay=1,
                            ann, 
                            maximalDropout,
                            gene.pair,
                            additivityMethod=c("bliss", "naive")
){
    additivityMethod <- match.arg(additivityMethod)
    x <- strsplit(gene.pair, split = "_")[[1]]
    gene1 <- x[[1]]
    gene2 <- x[[2]]
    sko1 <- .getSkoDataBX(se=ratios,
                          gene=gene1,
                          assay=assay,
                          sko.position="first")
    sko2 <- .getSkoDataBX(se=ratios,
                          gene=gene2,
                          assay=assay,
                          sko.position="second")
    sko1 <- mean(sko1$logRatios, na.rm=TRUE)
    sko2 <- mean(sko2$logRatios, na.rm=TRUE)
    if (additivityMethod=="naive"){
        dko_expected <- sko1+sko2
    } else if (additivityMethod=="bliss"){
        bliss_factor <- .getBlissFactor(sko1, 
                                       maximalDropout=maximalDropout)

        dko_expected <- sko1 + sko2*(bliss_factor)
    }
    return(dko_expected)
}

# Needed in Piru's dual guides analysis
.extractAllLFCs <- function(se,
                            gene.pair,
                            assay=1
){
    x <- strsplit(gene.pair, split = "_")[[1]]
    gene1 <- x[[1]]
    gene2 <- x[[2]]
    sko1 <- .getSkoDataBX(se=se,
                          gene=gene1,
                          assay=assay,
                          sko.position="first")
    sko2 <- .getSkoDataBX(se=se,
                          gene=gene2,
                          assay=assay,
                          sko.position="second")
    dko <- .getDkoDataBX(se=se,
                         gene1=gene1,
                         gene2=gene2,
                         assay=assay)
    y <- c(sko1$logRatios,
           sko2$logRatios,
           dko$logRatios)
    return(y)
}





# #' @title Generate Boxplot with scatterplot overlay of two given genes of interest from dualGuideDLFC logr output
# #' @description Generate Boxplot/scatterplot of two given genes of interest from dualGuideDLFC logr output.
# #' @param se.ratios Either a numeric matrix or a SummarizedExperiment.
# #' @param gene.pairs List of gene pairs.
# #' @param control Character string specifying "neg" or "NTC".
# #' @param xlim should be number of boxplots (guides) to be plotted. default is 11.
# #' @param ylim y axis limits, the function gathers this automatically.
# #' @param aspect.ratio default is 2/4, which is suitable for most plots.
# #' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_jitter position_jitter 
# #' @importFrom ggplot2 theme theme_bw element_text geom_vline geom_hline
# #' @export
# drawDualGuideBoxplotsPilot <- function(se,
#                                 gene1,
#                                 gene2,
#                                 assay=1,
#                                 xlim=NULL,
#                                 ylim=NULL,
#                                 aspect.ratio=2/5,
#                                 ...){
#   sko1     <- .getSkoDataBX(se=se, gene=gene1, assay=assay, sko.position="first", ...)
#   sko2     <- .getSkoDataBX(se=se, gene=gene2, assay=assay, sko.position="second", ...)
#     dko      <- .getDkoDataBX(se=se, gene1=gene1, gene2=gene2, assay=assay, ...)
#   all.data <- rbind(sko1, dko, sko2)
#   col      <- as.numeric(as.factor(all.data$col))
#   # Try either "Guide" or .data$Guide
#   # changed aes to aes_string
#   # vline1 <- length(sko1) + 0.5 # auto detect vline before target_target
#   # vline2 <- vline1 + 1 # auto detect vline after target_target
#   ggplot(all.data, aes_string(x="Guide", y="logRatios", Guide=FALSE)) +
#   geom_boxplot(outlier.shape=NA) +
#   geom_jitter(position=position_jitter(width=.2, height=0), color=col) +
#   theme_bw() +
#   theme(aspect.ratio=2/5, axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
#   geom_vline(xintercept=5.5) +
#   geom_vline(xintercept=6.5) +
#   geom_hline(yintercept=0)
# }


#' @title Generate within-sample normalization plots
#' @description Generate within-sample normalization plots.
#' @param se A SummarizedExperiment object.
#' @param sko.control String specifying the control non-cutting gene
#'     used in single-knockout constructs.
#' @param gene1.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 1.
#' @param gene2.field String specifying colum name in \code{rowData(se)} 
#'     containing gene name for gRNA in position 2.
#' @param gene.pair.field String specifying colum name in \code{rowData(se)}
#'     containing gene pair name.
#' @author Jean-Philippe Fortin
#' 
#' @export
drawWithinSampleNormalizationDensities <- function(se,
                                                   mode="all",
                                                   type.field="class",
                                                   type.levels.dko="or_or",
                                                   type.levels.sko1="or_neg",
                                                   type.levels.sko2="neg_or",
                                                   gene1.field="gene_symbol_1",
                                                   gene2.field="gene_symbol_2",
                                                   gene.pair.field="group",
                                                   sko.control="neg",
                                                   mfrow=c(3,8)
){
    assays(se)[[1]] <- as.matrix(assays(se)[[1]])
    Y <- assays(se)[[1]]
    Y <- log2(Y+1)
    if (mode=="all"){
        sko1.indices <- .getSkoIndices(se=se,
                                       sko.position="first",
                                       sko.control=sko.control,
                                       gene1.field=gene1.field,
                                       gene2.field=gene2.field,
                                       gene.pair.field=gene.pair.field)
        sko2.indices <- .getSkoIndices(se=se,
                                       sko.position="second",
                                       sko.control=sko.control,
                                       gene1.field=gene1.field,
                                       gene2.field=gene2.field,
                                       gene.pair.field=gene.pair.field)
        dko.indices <- .getDkoIndices(se=se,
                                      sko.control=sko.control,
                                      gene1.field=gene1.field,
                                      gene2.field=gene2.field,
                                      gene.pair.field=gene.pair.field)
    } else {
        dko.indices <- which(rowData(se)[[type.field]]==type.levels.dko)
        sko1.indices <- which(rowData(se)[[type.field]]==type.levels.sko1)
        sko2.indices <- which(rowData(se)[[type.field]]==type.levels.sko2)
    }
    par(mfrow=mfrow)
    for (k in 1:ncol(Y)){
        xlab="LogFCs"
        ylab="Density"
        densities <- list()
        densities[[1]] <- density(Y[sko1.indices,k], na.rm=TRUE)
        densities[[2]] <- density(Y[sko2.indices,k], na.rm=TRUE)
        densities[[3]] <- density(Y[dko.indices,k], na.rm=TRUE)
        plot(densities[[1]],
             col=2,
             main=colnames(Y)[k],
             xlab=xlab,
             ylab=ylab) #sko1
        lines(densities[[2]], col=3) #sko2
        lines(densities[[3]], col=1)
    }
}






