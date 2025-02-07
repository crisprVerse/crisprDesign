#' @title Draw relationship between coverage and a covariate.
#' @description Draw scatterplot of coverage vs a covariate.
#' @param object A SummarizedExperiment.
#' @param cutoff Numeric value specifying minimal number of reads required
#'     to call a barcode.
#' @param cov Numeric vector specifying covariate to visualize. Length must be
#'     equal to the number of columns in object.
#' @param xlab Character string of X axis label.
#' @param ref.field Character string specifying column in colData(se) 
#'     that contains the reference sample information.
#' @param ref.label Character string specifyng label of the reference sample
#'     to be found in colData(se)[[ref.field]].
#' @param pch Numeric value indicating the plot character for the points.
#' @param main Usual main argument in base plot.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @importFrom SummarizedExperiment assays colData
#' @export
plotCovariateRelationship <- function(object,
                                      cutoff=10,
                                      cov,
                                      ref.field="Group",
                                      ref.label=NULL,
                                      xlab="Covariate of interest",
                                      main="Title",
                                      pch=20,
                                      ...){
    key   <- .getKeyBarcoding(object, colname=ref.field)
    group <- colData(object)[[ref.field]]
    if (is.null(ref.label)){
        wh <- seq_len(ncol(object))
    } else {
        wh <- which(!group %in% ref.label)
    }
    key   <- key[wh,]
    Y     <- assays(object)[[1]][,wh]
    nbarcodes <- colSums(Y>=cutoff) 
    col <- key$col
    cov <- cov[wh,drop=FALSE]
    plot(cov, 
         nbarcodes,
         col=col, 
         pch=pch, 
         xlab=xlab,
         ylab=paste0("Number of unique barcodes with at least ",cutoff, " read(s)"),
         main=main,
         ...)
    .addScatterplotBarcodingLegend(key, ...)
}



#' @title Draw barplots of total number of retrieved barcodes. 
#' @description Draw barplots of total number of retrieved barcodes. 
#' @param object SummarizedExperiment object.
#' @param cutoff Minimum of reads required to detect a barcode.
#' @param sort Logical value, if TRUE samples will be drawn by color order.
#' @param col Optional character vector specifying sample color. Must be 
#'      of the same length as the number of columns in the provided object.
#' @param col.group.field Character string indicating the group field that
#'     includes the group label.
#' @param col.group.label Character string indicating the specific group to plot
#'      of the same length as the number of columns in the provided object.
#' @param excl.group.label Character string indicating the name of the group 
#'      you do not want to include in the plot (usually reference).
#' @param ylim Usual ylim argument in base plot.
#' @param xlab Usual xlab argument in base plot.
#' @param ylab Usual ylab argument in base plot.
#' @param main Usual main argument in base plot.
#' @param cex.names Numeric value indicating text size for labels.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @export
#' @importFrom SummarizedExperiment assays colData
drawBarplotsBarcoding <- function(object,
                                  cutoff=1,
                                  sort=TRUE,
                                  col=NULL,
                                  col.group.field=NULL,
                                  col.group.label=NULL,
                                  excl.group.label=NULL,
                                  ylim=c(0,10000),
                                  xlab="Sample",
                                  ylab="Number of barcodes",
                                  main="",
                                  cex.names=0.5,
                                  ...){
    col.null <- "black"
    group    <- colData(object)[[col.group.field]]
    if (!is.null(excl.group.label)){
        wh <- which((group %in% col.group.label) & !(group %in% excl.group.label))
    } else {
        wh <- which(group %in% col.group.label)
    }
    if (is.null(col)){
        pheno <- colData(object)
        if (!col.group.field %in% colnames(pheno)){
            stop("col.group.field is not found in colData(object)")
        }
        key <- .getKeyBarcoding(object, col.group.field)
        key <- key[wh,,drop=FALSE]   
        col <- key$col
    }
    if (is.null(col)) {
        col <- rep(col.null, ncol(counts))
    }
    object <- object[,wh,drop=FALSE]
    n      <- ncol(object)
    if(is(object, "SummarizedExperiment")){
        counts <- as.matrix(assays(object)[[1]])
    } else {
        stop("This function only supports 'SummarizedExperiment' objects.")
    }
    ys <- apply(counts, 2, function(x){
        x[x>=cutoff]
    })
    o <- rep(1, n)
    if (sort){
        o <- order(col)
    }
    ns <- vapply(ys, length, FUN.VALUE=0)
    barplot(ns[o],
            col=col[o],
            border=col[o],
            las=2,
            ylim=ylim,
            main=main,
            cex.names=cex.names,
            ...)
    .addBarplotsBarcodingLegend(key, ...)
}


#' @title Draw densities to represent barcode distributions.
#' @description Draw densities to represent barcode distributions.
#' @param object A SummarizedExperiment object.
#' @param log.scale Should the data be log-transformed? TRUE by default.
#' @param cutoff Minimum of reads required to detect a barcode.
#' @param normalize Should counts be normalized?
#' @param normalize_fun String specifying method (median or mean) to normalize counts.
#' @param col.group.field Character string indicating the field containing reference group.
#' @param col.group.label Character string indicating the name of the treatment group to color separately.
#' @param sort Should samples be drawn by color order? 
#' @param col Optional character vector specifying sample color. Must be 
#'      of the same length as the number of columns in the provided object.
#' @param lty Optional character vector specifying sample line type. Must be 
#'      of the same length as the number of columns in the provided object.
#' @param xlim Usual xlim argument in base plot.
#' @param ylim Usual ylim argument in base plot.
#' @param xlab Usual xlab argument in base plot.
#' @param ylab Usual ylab argument in base plot.
#' @param main Usual main argument in base plot.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @importFrom SummarizedExperiment assays colData
#' @importFrom matrixStats colMedians
#' @export
drawDensitiesBarcoding <- function(object,
                                   log.scale=TRUE,
                                   cutoff=1,
                                   normalize=TRUE,
                                   normalize_fun=c("mean", "median"),
                                   col.group.field=NULL,
                                   col.group.label=NULL,
                                   sort=TRUE,
                                   col=NULL, 
                                   lty=NULL,
                                   xlim=NULL,
                                   ylim=NULL,
                                   xlab="Log2 Barcode Counts", 
                                   ylab="Density",
                                   main="",
                                   ...){
    col.null      <- "black"
    col.fill      <- "white"
    normalize_fun <- match.arg(normalize_fun)
    if (is.null(lty)) {
        lty <- rep(1, ncol(object))
    }
    if (is(object,"SummarizedExperiment")){
        counts <- as.matrix(assays(object)[[1]])
    } else {
        stop("Only SummarizedExperiment objects are supported.")
    }
    if (is.null(col.group.label)) {
        wh  <- seq_len(ncol(object))
    } else {
        group <- colData(object)[[col.group.field]]
        wh    <- which(group %in% col.group.label)
    }
    counts <- counts[,wh,drop=FALSE]

    ### Normalization
    if (normalize){
        temp <- counts
        temp[temp==0] <- NA
        if (normalize_fun=="median"){
            factors <- colMedians(temp, na.rm=TRUE) 
        } else {
            factors <- colMeans(temp, na.rm=TRUE)   
        }
        factors <- factors/median(factors, na.rm=TRUE)
        counts <- sweep(counts,2,factors, "/")
    }
    ys <- apply(counts,2, function(x){
        x[x>=cutoff,drop=FALSE]
    })
    ### Log-scale
    if (log.scale){
        ys <- lapply(ys, function(x) log2(x+1))
    }
    n <- ncol(counts)
    densities <- lapply(seq_len(n), function(i){
        density(ys[[i]], na.rm=TRUE)
    })
    minx   <- min(unlist(lapply(densities, function(x) min(x$x))))
    maxx   <- max(unlist(lapply(densities, function(x) max(x$x))))
    miny   <- min(unlist(lapply(densities, function(x) min(x$y))))
    maxy   <- max(unlist(lapply(densities, function(x) max(x$y))))
    rangex <- abs(maxx - minx)
    rangey <- abs(maxy - miny)
    epsx   <- 0.1 * rangex
    epsy   <- 0.1 * rangey
    if (is.null(xlim)) {
        xlim = c(minx - epsx, maxx + epsx)
    }
    if (is.null(ylim)) {
        ylim = c(miny, maxy + epsy)
    }
    if (is.null(col)){
        pheno <- colData(object)
        if (!col.group.field %in% colnames(pheno)){
            stop("col.group.field is not found in colData(object)")
        }
        key <- .getKeyBarcoding(object, col.group.field)
    } 
    key <- key[wh,,drop=FALSE]
    col <- key$col
    o <- rep(1, n)
    if (sort){
        o <- order(col)
    }
    plot(densities[[1]], 
         col = col.fill, 
         bty = "L", 
         xlim = xlim, 
         ylim = ylim, 
         xlab = xlab, 
         ylab = ylab, 
         main = main, 
         ...)
    for (i in seq_len(ncol(counts))) {
        lines(densities[o][[i]], 
              col = col[o][i], 
              lty = lty[o][i], 
              ...)
    }
    .addDensitiesBarcodingLegend(key,...)
}







#' @title Draw number of reads and number of unique barcodes. 
#' @description Draw a scatterplot showing number of reads and number of unique barcodes given a coverage cutoff. 
#' @param object A Summarized Experiment.
#' @param col.group.field Character string indicating the field containing reference group.
#' @param col.group.label Character string indicating the name of the treatment group to color separately.
#' @param excl.group.label Character string indicating the name of the group you do not want to include in the plot (usually reference).
#' @param cutoff Numeric value.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @export
plotCoverageRelationship <- function(object,
                                     col.group.label=NULL,
                                     col.group.field=NULL,
                                     excl.group.label=NULL,
                                     cutoff=1,
                                     ...){
    xlab  <- "Number of reads"
    pch   <- 20
    ylab  <- paste0("Number of unique barcodes with at least ",cutoff, " read(s)")
    group <- colData(object)[[col.group.field]]
    if (!is.null(excl.group.label) & !is.null(col.group.label)){
        wh <- which((group %in% col.group.label) & !(group %in% excl.group.label))
    } 
    if (!is.null(excl.group.label) & is.null(col.group.label)){
        wh <- which(!group %in% excl.group.label)
    } else {
        wh <- which(group %in% col.group.label)
    }
    pheno <- colData(object)
    Y  <- assays(object)[[1]][,wh,drop=FALSE]
    if (!"nreads" %in% names(pheno)){
        stop("'nreads' not found as a column in colData")
    }
    nreads <- colSums(Y)
    nbarcodes <- colSums(Y>=cutoff)
    key     <- .getKeyBarcoding(object, col.group.field)
    key     <- key[wh,,drop=FALSE]
    col     <- key$col
    plot(nreads,
         nbarcodes,
         col=col,
         pch=pch,
         xlab=xlab,
         ylab=paste0("Number of unique barcodes with at least ",cutoff, " read(s)"),
         ...)
    .addScatterplotBarcodingLegend(key,...)
}







#' @title Draw distributions of barcode clone sizes and number of cells per barcode.
#' @description Draw distributions of barcode clone sizes and number of cells per barcode.
#' @param y Numeric vector of barcode counts.
#' @param ylim Usual ylim argument in base plot.
#' @param start Numeric value to specifyc start bin on log2 scale.
#' @param end Numeric value to specifyc start bin on log2 scale.
#' @param step Numeric value to specify step value to control number of bars.
#' @param normFactor Numeric value to control number of reads per clone 
#' @param main Usual main argument in base plot.
#' @param cex Numeric value indicating size of labels.
#' @param cex.axis Numeric value indicating size of axis labels.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @export
drawDualDistributionsBarcoding <- function(y,
                                           start=4,
                                           end=25,
                                           step=0.1,
                                           normFactor=2000,
                                           ylim=c(0,300000),
                                           main="",
                                           cex=0.7,
                                           cex.axis=0.4,
                                           ...){
    col   <- c("deepskyblue3", "orange")
    cats1 <- seq(start,end,step)
    cats1 <- cats1[-length(cats1)]
    cats2 <- c(cats1[-1],end)
    bins  <- cbind(cats1,cats2)
    nbarcodes <- lapply(seq_len(nrow(bins)), function(i){
        sum(y>=bins[i,1] & y<bins[i,2])
    })
    nreads <- lapply(seq_len(nrow(bins)), function(i){
        wh <- which(y>=bins[i,1] & y<bins[i,2])
        sum(2^y[wh]-1)
    })
    nreads     <- unlist(nreads)
    nbarcodes  <- unlist(nbarcodes)
    nreads.cum <- cumsum(nreads)
    nreads.cum <- nreads.cum/sum(nreads)
    data       <- cbind(nbarcodes, nreads/normFactor)
    bins2      <- 2^bins
    meds       <- (bins2[,2]+bins2[,1])/2
    meds       <- round(meds)
    barplot(t(data), 
            beside=TRUE, 
            main=main,
            col=col,
            xlab="Colony size (number of reads per barcode)",
            ylab="",
            border=NA,
            ylim=ylim,
            ...)
    results <- barplot(t(data),
                       beside=TRUE,
                       plot=FALSE,
                       ...)
    x <- (results[1,]+results[2,])/2
    axis(side=1, 
         at=x, 
         meds, 
         las=2, 
         cex.axis=cex.axis,
         ...)
    legend("topleft", 
           fill=col, 
           border=col, 
           c("Number of unique barcodes", 
           paste0("Number of reads divided by ", normFactor)),
           cex=cex, 
           bty="n",
           ...)
}




#' @title Draw cumulative distributions of barcode clone sizes and number of cells per barcode.
#' @description Draw cumulative distributions of barcode clone sizes and number of cells per barcode.
#' @param y Numeric vector of barcode counts.
#' @param ylim Usual ylim argument in base plot.
#' @param start Numeric value to specifyc start bin on log2 scale.
#' @param end Numeric value to specifyc start bin on log2 scale.
#' @param step Numeric value to specify step value to control number of bars.
#' @param normFactor Numeric value to control number of reads per clone.
#' @param main Usual main argument in base plot.
#' @param cex Numeric value indicating size of labels.
#' @param cex.axis Numeric value indicating size of axis labels.
#' @param ... Other arguments to passe to base plot function.
#' @return Returns nothing. A plot is produced as a side effect. 
#' @importFrom graphics axis
#' @export
drawCumulativeDistributionsBarcoding <- function(y,
                                                start=6, 
                                                end=20, 
                                                step=0.1, 
                                                normFactor=2000,
                                                ylim=c(0,250000), 
                                                main="",
                                                cex.axis=0.4,
                                                cex=0.7,
                                                ...){
    col       <- c("deepskyblue3", "orange")
    cats1     <- seq(start,end,step)
    cats1     <- cats1[-length(cats1)]
    cats2     <- c(cats1[-1],end)
    bins      <- cbind(cats1,cats2)
    nbarcodes <- lapply(seq_len(nrow(bins)), function(i){
        sum(y>=bins[i,1] & y<bins[i,2])
    })
    nreads <- lapply(seq_len(nrow(bins)), function(i){
        wh <- which(y>=bins[i,1] & y<bins[i,2])
        sum(2^y[wh]-1)
    })
    nreads    <- unlist(nreads)
    nbarcodes <- unlist(nbarcodes)
    nbarcodes <- cumsum(nbarcodes)/sum(nbarcodes)
    nreads    <- cumsum(nreads)/sum(nreads)
    data      <- cbind(nbarcodes, nreads)
    bins2     <- 2^bins
    meds      <- (bins2[,2]+bins2[,1])/2
    meds      <- round(meds)
    barplot(t(data), 
            beside=TRUE, 
            main=main,
            col=col,
            xlab="Colony size (number of reads per barcode)",
            ylab="",
            border=NA, 
            ylim=ylim,
            ...)
    results <- barplot(t(data), 
                       beside=TRUE, 
                       plot=FALSE)
    x <- (results[1,]+results[2,])/2
    axis(side=1, 
         at=x, 
         meds, 
         las=2, 
         cex.axis=cex.axis,
         ...)
    legend("topleft", 
           fill=col, 
           border=col, 
           c("Number of unique barcodes", paste0("Number of reads divided by ", normFactor)),
           cex=cex, 
           bty="n",
           ...)
}




# # Looking at densities for reference and expansion:
# addDensitiesBarcodingLegend <- function(key, where="topright",cex=0.75){
#   key <- key[!duplicated(key),]
#   legend(where, col=key$col, legend=key$label,
#      lty=1, bty='n',cex=cex
#   )
# }

# addBarplotsBarcodingLegend <- function(key, where="topright",cex=0.75){
#     key <- key[!duplicated(key),]
#     legend(where, col=key$col, fill=key$col, border=key$col,
#      legend=key$label, bty='n',cex=cex
#     )
# }



#' @importFrom SummarizedExperiment colData
.getKeyBarcoding <- function(object, 
                             colname="Group"){
    names <- names(colData(object))
    if (!colname %in% names){
        stop("colname column cannot be found in colData(se)")
    }
    group <- colData(object)[[colname]]
    key <- data.frame(col=as.numeric(as.factor(group)),
                      stringsAsFactors=FALSE)
    key$label <- group
    return(key)
}



.addScatterplotBarcodingLegend <- function(key, 
                                           where="topleft", 
                                           cex=0.75,
                                           pch=20,
                                           bty="n",
                                           ...){
    key <- key[!duplicated(key),,drop=FALSE]
    legend(where, 
           col=key$col, 
           pch=pch,
           legend=key$label, 
           bty=bty,
           cex=cex,
           ...)
}


.addDensitiesBarcodingLegend <- function(key, 
                                         where="topright",
                                         cex=0.75,
                                         bty="n",
                                         ...){
    key <- key[!duplicated(key),,drop=FALSE]
    legend(where, 
           col=key$col, 
           legend=key$label,
           lty=1, 
           bty=bty,
           cex=cex,
           ...)
}


.addBarplotsBarcodingLegend <- function(key, 
                                        where="topright",
                                        cex=0.75,
                                        bty='n',
                                        ...){
    key <- key[!duplicated(key),,drop=FALSE]
    legend(where, 
           col=key$col, 
           fill=key$col, 
           border=key$col,
           legend=key$label, 
           bty=bty,
           cex=cex,
           ...)
}
