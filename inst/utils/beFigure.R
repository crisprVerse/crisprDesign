getGenes <- function(x){
    strsplit(x, split="_") %>% lapply(., function(x) x[[1]]) %>% unlist
}


normalizeNTC <- function(Z){
  factors <- colMedians(Z[grepl("NTC", rownames(Z)),])
  factors <- factors-median(factors)
  Z <- sweep(Z, 2, factors, "-")
  a <- min(Z)
  if (a<0){
    Z <- Z+abs(a)
  }
  Z
}

normalizeTotal <- function(Z){
  factors <- colMeans(Z)
  factors <- factors-median(factors)
  Z <- sweep(Z, 2, factors, "-")
  Z
}


# To get annotation;
library(crisprDesign)
library(crisprDesignAux)
library(crisprDesignGne)
library(crisprDesignData)
bsgenome <- getGenomePackage()
library(SummarizedExperiment)
library(dplyr)
load("../objects/resultsBase.rda")
se <- resultsBase[["se"]]
options(tibble.width = Inf)

# Getting library annotation:
ann <- rowData(se) %>% as.data.frame
ann <- dplyr::rename(ann, guide_19mer=barcode)
ann <- dplyr::rename(ann, id=ID)
anns <- split(ann, f=ann$gene_symbol)
genes <- c("myGene")
anns <- anns[genes]
cut.offset=-15 #Optimal window around -15
nuc <- SpCas9
for (i in 1:6){
    nuc@cutSites[i] <- cut.offset    
}

#Generating alignments and adding tx info:
anns <- lapply(anns, function(ann){
    ann$cut_site <- crisprBase:::getCutSiteFromPamSite(strand=ann$strand, 
        pam_site=ann$pam_site, 
        nuclease=nuc
    )
    ann
})
ann <- anns[[1]]
ann$pam <- "AGG"
ann <- as.data.frame(ann)

gs <- crisprDesignAux::DataFrameToGuideSet(df=ann,
                                           bsgenome=bsgenome,
                                           crisprNuclease=SpCas9,
                                           spacerCol="guide_19mer",
                                           idCol="id")
gs <- addGeneAnnotation(gs, txObject=txdb_human)
geneAnn <- geneAnnotation(gs)
txId <- "myTxId"
geneAnn <- geneAnn[geneAnn$tx_id==txId,]
stopifnot(!any(duplicated(rownames(geneAnn))))
gs$aa <- NA
wh <- match(rownames(geneAnn), names(gs))
gs$aa[wh] <- geneAnn$aminoAcidIndex
gs$geneAnnotation <- NULL
out <- GuideSet2DataFrames(gs)[[1]]
mainProteinAnnotation <- out
save(mainProteinAnnotation,
     file="../objects/mainProteinAnnotation.rda")



library(crisprDesign)
library(crisprScreen)
library(stringr)
library(dplyr)
library(Biobase)
library(SummarizedExperiment)
library(matrixStats)
library(RColorBrewer)
library(scales)
library(graphics)
library(limma)
library(plotly)
library(magrittr)
filteringCutoff <- 100
logFCCutoff <- 1
includeDay0 <- TRUE
source("utils.R")




load("../objects/resultsBase.rda")
se <- resultsBase[["se"]]
stats <- resultsBase[["barcode"]]



getDomains <- function(){
    temp <- read.table(paste0("../data/domains_","gene","_pfam.txt"))
    colnames(temp) <- c("source", "name", "start", "end")
    temp <- temp[temp$source=="Pfam",]
    temp$col <- 1:nrow(temp)
    domains <- temp
    domains$start %<>% as.integer
    domains$end %<>% as.integer
    return(domains)
}
domains <- getDomains()



load("../objects/mainProteinAnnotation.rda")
mainProteinAnnotation <- as.data.frame(mainProteinAnnotation)
aaAnn <- mainProteinAnnotation
dfs <- list(aaAnn)
dfs <- lapply(dfs, function(x){
    rownames(x) <- x$ID
    x
})



dfs <- lapply(dfs, function(x){
    x[!is.na(x$aa),]
})




# Let's add domains to the annotations:
dfs <- lapply(1:length(dfs), function(i){
    df <- dfs[[i]]
    domains$start %<>% as.integer
    domains$end %<>% as.integer
    gr.domains <- GRanges("chr1", IRanges(start=domains$start, end=domains$end))
    
    gr.guides  <- GRanges("chr1", IRanges(start=df$aa, end=df$aa))
    wh <- findOverlaps(gr.guides, gr.domains)
    wh1 <- queryHits(wh)
    wh2 <- subjectHits(wh)
    df$domain <- NA
    df$col <- "grey75"
    df$domain[wh1] <- domains[wh2,"name"]
    df$col[wh1] <- domains[wh2,"col"]
    dfs[[i]] <- df
})




proteinAnn <- dfs[[1]]


pToAlpha <- function(p,
                     cats1=c(0.005, 0.01, 0.05)
){
    color <- p
    color[p > cats1[3]] <- 0.1
    color[p <= cats1[3]] <- 1
    color[p <= cats1[2]] <- 1
    color[p <= cats1[1]] <- 1
    color
}


haleyLabPalette <- function(){
    codes <- c("#001219",
               "#005F73",
               "#0A9396",
               "#94D2BD",
               "#E9D8A6",
               "#EE9B00",
               "#CA6702",
               "#BB3E03",
               "#AE2012",
               "#9B2226")
    codes
}



coloring <- "replicate"
ylim=c(-1,10)

plotLFC <- function(compName.list,
                    topCond.list,
                    bottomCond.list,
                    ylab="logFC",
                    main="",
                    ylim=c(-1,10),
                    xlim=NULL,
                    coloring=c("replicate","condition"),
                    saveLFCs=FALSE,
                    legend=TRUE
){

    coloring <- match.arg(coloring)
    ann <- proteinAnn    
    xlab <- ""
    ylab="Log2 fold change"

    if (coloring=="replicate"){
        col1 <- haleyLabPalette()[4]
        col2 <- haleyLabPalette()[2]
        col3 <- haleyLabPalette()[6]
        col <- c(col1,col2,col3)
        #col <- ifelse(grepl("RepA", samples2), colA, colB)
    } else {
        col1 <- haleyLabPalette()[2]
        col2 <- haleyLabPalette()[6]
        col <- c(col1, col2)
    }



    for (kk in 1:length(compName.list)){
        compName <- compName.list[[kk]]
        topCond <- topCond.list[[kk]]
        bottomCond <- bottomCond.list[[kk]]
        pheno <- colData(se)
        pval <- stats[[compName]]$PValue
        beta <- stats[[compName]]$LogFC
        fdr  <- stats[[compName]]$FDR
        names(pval) <- names(beta) <- names(fdr) <- rownames(stats[[compName]])
        pval[beta<=0] <- 1
        fdr[beta<=0]  <- 1

        Z <- assays(se)[[2]] #logcounts
        y <- Z[, pheno$Condition==topCond]
        x <- Z[, pheno$Condition==bottomCond]
        lfcs <- y-x

        lfcs <- lfcs[rownames(ann),]
        beta <- beta[rownames(ann)] 
        fdr <- fdr[rownames(ann)]
        pval <- pval[rownames(ann)]

        #ALpha
        #pp <- pval[, grepl(drug, colnames(pval))]
        pp <- fdr
        alpha <- pToAlpha(pp)
        alpha <- alpha[rownames(ann)]

        x <- ann$aa
        Y <- lfcs
        B <- beta

        if (kk==1){
            plot(x, Y[,1],main=main, xlab=xlab, ylab=ylab, col="white", ylim=ylim, xlim=xlim)
        }

        if (grepl("tx1", topCond)){
            col <- col1
        } else {
            col <- col2
        }
        for (kkk in 1:ncol(Y)){
            if (coloring=="replicate"){
                current.col <- rep(col[kkk], length(x))
            } else {
                current.col <- rep(col, length(x))
            }
            
            current.col[which(B<=logFCCutoff)] <- "grey75"
            current.col[which(alpha==0.1)] <- "grey75"
            current.col[is.na(alpha)] <- "white"
            current.col <- alpha(current.col, alpha)
            current.col[is.na(alpha)] <- alpha(current.col[is.na(alpha)],0)
            current.col[which(B<=logFCCutoff)] <- alpha(current.col[which(B<=logFCCutoff)],0.1)
            points(x,
                   Y[,kkk],
                   col=current.col,
                   pch=20,
                   cex=1.2)
        }
    }

    abline(h=0, lty=3)
    if (legend){
        if (coloring=="replicate"){
            legend("topleft",
                   col=c(col1, col2, col3), 
                   c("Rep1", "Rep2", "Rep3"), pch=20)
        } else {
            legend("topleft", 
                   col=c(col1,col2),
                   c("tx2", "tx2"), pch=20)
        }
    }
}



addDomains <- function(){

    plotDomains <- function(domains, ypos, w=0.1){
        for (i in 1:nrow(domains)){
          domain <- domains[i,]
          rect(xleft=domain[,"start"],
               xright=domain[,"end"],
               ybottom=ypos-w,
               ytop= ypos+w,
               col = domains[i, "col"],
               border = domains[i,"col"]
          )
        }
    }

    addDomainNames <- function(domains, ypos, w=0.1, cex=1){
        for (i in 1:nrow(domains)){
            label <- domains[i, "name"]
            mid <- (domains[i,"start"]+domains[i,"end"])/2
            text(x=mid, y=ypos, label, cex=cex)
        }
    }
    domainsi <- domains
    newcol <- rep(c("grey55", "grey75"), nrow(domainsi))
    newcol <- newcol[1:nrow(domainsi)]
    domainsi$col <- newcol
    dfi <- proteinAnn
    xx <- 0:(max(dfi$aa)+10)
    plot(xx,
         rep(1,length(xx)),
         bty="n",
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         col="white",
         ylim=c(0,2))

    y.domains <- 1
    segments(x0=0,
             x1=max(dfi$aa),
             y0=y.domains,
             y1=y.domains)
    plotDomains(domainsi,
                ypos=y.domains,
                w=0.3)
    addDomainNames(domainsi,
                   ypos=0.2,
                   w=0.5,
                   cex=0.7)
}




createFinalFigure <- function(compName.list,
                              topCond.list,
                              bottomCond.list,
                              coloring=c("replicate","condition"),
                              ylim,
                              legend=TRUE,
                              main=""
){
    coloring <- match.arg(coloring)
    par(fig=c(0,1,0.15,1), new=FALSE)
    plotLFC(compName.list=compName.list,
            topCond.list=topCond.list,
            bottomCond.list=bottomCond.list,
            coloring=coloring,
            ylim=ylim,
            legend=legend,
            main=main)
    grid()
  
    par(fig=c(0,1,0,0.4), new=TRUE)
    addDomains()
    mtext(text = "Amino Acid Position",
          side = 1,#side 1 = bottom
          line = 1)
}












