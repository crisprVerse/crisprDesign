#' @title Generate plot of barcode distributions on a log2 scale (static plot)
#' @description Generate plot of barcode distributions on a log2 scale (static plot).
#' @param object Either a numeric matric or a SummarizedExperiment
#' @param log.scale Should the data be log-transformed? TRUE by default.
#' @param sort Should samples be drawn by color order? 
#' @param col Optional character vector specifying sample color. Must be 
#' @param col.group.field Optional string indicating which column of colData(object) should be
#' used to color densities. 
#' @param lty Optional character vector specifying sample line type. Must be 
#' of the same length as the number of columns in the provided object.
#' @param xlim Usual xlim argument in base plot
#' @param ylim Usual ylim argument in base plot
#' @param xlab Usual xlab argument in base plot
#' @param ylab Usual ylab argument in base plot
#' @param main Usual main argument in base plot
#' @param legend.where Where should be the legend?
#' @param assay_index Integer specifying which assay should be used. 1 by default. 
#' @param ... Other arguments to passe to base plot function
#' @return Returns nothing. A plot is produced as a side effect. 
#' @examples 
#' \dontrun{
#'  data("seExample")
#'  drawDensities(seExample)
#' }
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom matrixStats rowMedians
#' @importFrom stats density
#' @importFrom graphics lines plot
drawDensities <- function(object,
    log.scale=TRUE, 
    sort=TRUE, 
    col=NULL,
    col.group.field="Group",
    lty=NULL,
    xlim=NULL, 
    ylim=NULL, 
    xlab="Log2 sgRNA Counts", 
    ylab="Density", 
    main="",
    legend.where="topleft",
    assay_index=1,
    ...
){
    if (is.null(lty)){
    	lty <- rep(1, ncol(object))
    }
    n <- ncol(object)
	#if (class(object)=="ExpressionSet") counts <- exprs(object)
    if (is(object, "SummarizedExperiment")){
        xs <- assays(object)
        if (assay_index>length(xs)){
            stop("Assay index is greated than the number of assays. ")
        } else {
            counts <- as.matrix(xs[[assay_index]])
        }
    } else {
        counts <- object
    }
	if (log.scale){
        counts <- log2(counts+1)
    } 
	#if (is.null(col)){
    #    col <- rep("black", ncol(counts))
    #}
    if (is.null(col)){
        pheno <- colData(object)
        if (!col.group.field %in% colnames(pheno)){
            stop("col.group.field is not found in colData(object)")
        }
        key     <- data.frame(label=pheno[[col.group.field]])
        key$col <- as.numeric(as.factor(key$label))
        col <- key$col
    } 
	o <- rep(1, n)
	if (sort){
        o <- order(col)    
    }
    densities <- lapply(1:n, function(i){
        density(counts[,o][,i], na.rm=TRUE)
    })

    # Setting up axes limits:
    minx <- min(unlist(lapply(densities, function(x) min(x$x))))
    maxx <- max(unlist(lapply(densities, function(x) max(x$x))))
    miny <- min(unlist(lapply(densities, function(x) min(x$y))))
    maxy <- max(unlist(lapply(densities, function(x) max(x$y))))
    rangex <- abs(maxx-minx)
    rangey <- abs(maxy-miny)
    epsx <- 0.1*rangex
    epsy <- 0.1*rangey

    if (is.null(xlim)){
        xlim=c(minx-epsx,maxx+epsx)
    } 
    if (is.null(ylim)){
      ylim=c(miny,maxy+epsy)  
    } 

	plot(densities[[1]], col="white", bty="L", 
        xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main,...)
	for (i in 1:ncol(counts)){
		lines(densities[[i]], col=col[o][i], lty=lty[o][i],...)
	}
    legend(legend.where, col=unique(key$col), legend=unique(key$label), lty=1, bty="n")
}

#load("../objects/se.rda")
#se.norm <- normalizeNtc(se)
#col <- colData(se.norm)$Condition %>% as.factor %>% as.numeric
#drawDensities(se,col=col)
#drawDensities(se.norm,col=col)
#dev.off()

               



#' @title Compositional plot for ORF/Barcoding experiments.
#' @description Generate compositional barplots of ORF/barcode distributions. 
#' @param se Either a numeric matric or a SummarizedExperiment.
#' @param howmany How many ORFs should be shown and labeled on the plot?
#' @return Returns nothing. A plot is produced as a side effect. 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices palette
#' @importFrom utils capture.output
#' @importFrom graphics barplot legend lines par plot plot.new abline mtext
compositionalPlot <- function(se, 
                              howmany=10){
    palette(brewer.pal(howmany+1,"RdYlBu"))
    Y <- as.matrix(assays(se)[[1]])
    Y <- sweep(Y,2,colSums(Y),"/")
    A <- rowMeans(Y[,1:(ncol(Y)-1)])
    o <- order(A)
    Y <- Y[o,]

    n <- nrow(Y)
    Y <- Y[n:1,]
    Y <- rbind(Y[1:howmany,], colSums(Y[(howmany+1):n,]))
    if (!"Sample" %in% colnames(colData(se))){
        stop("`Sample` must be provided in colData(se)")
    }
    colnames(Y) <- colData(se)$Sample
    par(mar=c(10,4,4,4))
    layout(matrix(c(1,2), nrow = 1), widths = c(0.7, 0.5))
    barplot(Y*100, col=1:(howmany+1), ylab="Percentage of reads",las=2)
    plot.new()
    names <- rownames(Y)
    names[howmany+1] <- "Others"
    legend("topleft",fill=1:(howmany+1), names, bty="n", cex=1.1)
}






#load("../../ngsProjects/ngs3565_crispra_necroptosis/objects/resultsBase.rda")
#results <- resultsBase
#comp <- "Difference between `LPS` vs `NoLPS`"
#method="fry"
#par(mar=c(10,4,4,4))
#plotTopTargets(results, comp=comp, enrich=TRUE, method="rho", target=30)

#' @title Plot top targets.
#' @description Plots top targets.
#' @param results A results voom-screen object.
#' @param targets Numeric value specifying numer of targets to plot.
#' @param method Character string of results method (e.g. "fry", "simes", "holm-min", "rho").
#' @param enrich Logical value.
#' @param comp Character string specifying the comparison of interest.
#' @return Returns nothing. a plot is produced as a side effect.
#' @export
#' @importFrom SummarizedExperiment rowData
plotTopTargets <- function(results,
    targets=10,
    method="fry",
    enrich=TRUE,
    comp=NULL
){
    current.graphic.params <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(par(current.graphic.params)))

    if(!(enrich %in% c(TRUE, FALSE))){
        stop('enrich must be either TRUE or FALSE.')
    }
    
    if (!method %in% names(results[["gene"]])){
        stop(paste0("method ", method," not available."))
    }

    #Getting annotation
    ann <- rowData(results[["se"]])


    df      <- as.data.frame(results[["barcode"]][[comp]])
    df_gene <- as.data.frame(results[["gene"]][[method]][[comp]])
    wh <- match(df$group, df_gene$id)
    df$PValue_group <- df_gene$PValue[wh]
    df$Direction    <- df_gene$Direction[wh]
    df <- df[!is.na(df$PValue_group),]
    df <- df[!df$Direction=="mixed",]

    # Sorting
    if (enrich){
        df <- df[df$Direction=="up",,drop=FALSE]
        df <- df[order(df$PValue_group, df$group, -df$LogFC),]
        main <- "Enriched Targets"
    } else {
        df <- df[df$Direction=="down",,drop=FALSE]
        df <- df[order(df$PValue_group, df$group, df$LogFC),]
        main <- "Depleted Targets"
    }

    # Special treatment for Rho:
    if (method=="rho" & enrich){
        wh <- match(df$group, df_gene$id)
        df$Rho_group <- df_gene$Rho.enrichment[wh]
        df <- df[order(df$PValue_group,
            df$Rho_group, df$group, -df$LogFC),]
    }
    if (method=="rho" & !enrich){
        wh <- match(df$group, df_gene$id)
        df$Rho_group <- df_gene$Rho.depletion[wh]
        df <- df[order(df$PValue_group,
            df$Rho_group, df$group, df$LogFC),]
    }
   

    if(is.character(targets)){
      toptargets <- intersect(targets, ann$group)
      ntargets   <- length(toptargets)
      main  <- ""
    } else {
      if((length(targets) != 1) | !is.numeric(targets)){
        stop('"targets" must be specified as a single number or a vector of elements in the geneSymbol column of the annotation object.')
      }
      ntargets  <- targets
      main <- paste('Top', ntargets, main)
      toptargets <- unique(df$group)[1:ntargets]
    } 
    if(ntargets <= 0){
        stop("No valid targets were specified.")
    }
    
    targetrows <- which(df$group %in% toptargets)
    nguides <- unlist(lapply(toptargets, function(x){
        sum(df$group %in% x, na.rm = TRUE)
    }))

    lfc <- df$LogFC[targetrows]
    
    
    
    ylim <- c(min(0, min(lfc, na.rm=TRUE)), max(0, max(lfc, na.rm=TRUE)))
    
    #Compose a vector of the x locations
    x <- as.vector(unlist(mapply(function(start, numguides){if(numguides > 1){seq(start, start+1, length.out = numguides)} else{start + 0.5}}, 
                   start = seq(1, by = 2, length.out = ntargets), 
                   numguides = nguides, 
                   SIMPLIFY = TRUE)))
  
    
    
    #make the plot
    plot(x, lfc,
        xaxt='n',main = main, ylim = ylim,
        xlab = "",
        ylab = "Log2 gRNA Abundance Change",
        pch = 18,
        col = "firebrick2"
    )
    
    #points(xloc, lfc, pch = 18, col = "darkred")
    suppress <- capture.output(lapply(seq(2.5, by = 2, length.out = (ntargets - 1)), function(x){abline(v = x, lty = "dotted", col= "gray")}))
    abline(h = 0)
    #Add the labels
    axis(1, at = seq(1.5, by = 2, length.out = ntargets), labels = toptargets, las = 3)
    if (!is.null(comp)) {
        mtext(comp, side=3, line=0)
    }
}










