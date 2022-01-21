# #' @export
# getCbeWeights <- function(len=5, power=3){
#     if (len==5){
#         ws=c(1,1,1,2,3,4,3,2,1,1,1)
#         ws=ws^power/sum(ws^power)
#         names(ws) <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
#     }
#     if (len==4){
#         ws=c(1,1,2,3,4,3,2,1,1)
#         ws=ws^power/sum(ws^power)
#         names(ws) <- c(-4,-3,-2,-1,0,1,2,3,4)
#         ws <- ws
#     }
#     if (len==3){
#         ws=c(1,2,3,4,3,2,1)
#         ws=ws^power/sum(ws^power)
#         names(ws) <- c(-3,-2,-1,0,1,2,3)
#         ws <- ws
#     }
#     return(ws)
# }

# #' @export
# addCbeStopSites <- function(txAnn, cbe_type=c("c2t")){
#     cbe_type    <- match.arg(cbe_type)
#     gene_strand <- txAnn$strand[[1]]
#     stops  <- c("TAG","TAA","TGA")

#     temp <- txAnn[!is.na(txAnn$aa_id),,drop=FALSE]
#     temp$stop_fwd <- FALSE
#     temp$stop_rev <- FALSE
#     codons <- split(temp, f=temp$aa_id)
    
#     # On the forward strand:
#     codons <- lapply(codons, function(x){
#         y <- x$nuc
#         if (cbe_type=="c2t"){
#             y[y=="C"] <- "T"
#         }
#         codon <- paste0(y,collapse="")
#         if (codon %in% stops){
#             if (gene_strand=="+"){
#                 x$stop_fwd <- TRUE
#             } else {
#                 x$stop_rev <- TRUE
#             }
#         } 
#         x
#     })
#     codons <- lapply(codons, function(x){
#         y <- x$nuc
#         if (cbe_type=="c2t"){
#             y[y=="G"] <- "A"
#         }
#         codon <- paste0(y,collapse="")
#         if (codon %in% stops){
#             if (gene_strand=="+"){
#                 x$stop_rev <- TRUE
#             } else {
#                 x$stop_fwd <- TRUE
#             }
#         } 
#         x
#     })
#     df <- do.call(rbind, codons)
#     txAnn$stop_fwd <- txAnn$pos %in% df$pos[df$stop_fwd]
#     txAnn$stop_rev <- txAnn$pos %in% df$pos[df$stop_rev]
#     return(txAnn)
# }



# #' @importFrom Biostrings DNAString translate
# #' @export
# addEditedAlleles <- function(txAnn, cbe_type = c("c2t")){
#     cbe_type <- match.arg(cbe_type)
#     gene_strand <- txAnn$strand[[1]]
#     stops <- c("TAG", "TAA", "TGA")
#     temp <- txAnn[!is.na(txAnn$aa_id), , drop = FALSE]
#     dna <- paste0(temp$nuc, collapse="")
#     dna <- DNAString(dna)
#     if (cbe_type=="c2t"){
#         dna_edited_1 <- DNAString(gsub("C", "T", dna))
#         dna_edited_2 <- DNAString(gsub("G", "A", dna))
#     }
#     protein <- translate(dna)
#     protein_edited_1 <- translate(dna_edited_1)
#     protein_edited_2 <- translate(dna_edited_2)
#     protein <- strsplit(as.character(protein), "")[[1]]
#     protein_edited_1 <- strsplit(as.character(protein_edited_1), "")[[1]]
#     protein_edited_2 <- strsplit(as.character(protein_edited_2), "")[[1]]
#     protein <- rep(protein, each=3)
#     protein_edited_1 <- rep(protein_edited_1, each=3)
#     protein_edited_2 <- rep(protein_edited_2, each=3)
#     dna_edited_1 <- strsplit(as.character(dna_edited_1), "")[[1]]
#     dna_edited_2 <- strsplit(as.character(dna_edited_2), "")[[1]]
#     if (gene_strand=="+"){
#         temp$aa_edited_fwd  <- protein_edited_1
#         temp$aa_edited_rev  <- protein_edited_2
#         temp$nuc_edited_fwd <- dna_edited_1
#         temp$nuc_edited_rev <- dna_edited_2
#     } else {
#         temp$aa_edited_fwd  <- protein_edited_2
#         temp$aa_edited_rev  <- protein_edited_1
#         temp$nuc_edited_fwd <- dna_edited_2
#         temp$nuc_edited_rev <- dna_edited_1
#     }
#     temp$class_fwd <- NA
#     temp$class_fwd[temp$aa_edited_fwd==temp$aa] <- "silent"
#     temp$class_fwd[temp$aa_edited_fwd!=temp$aa] <- "missense"
#     temp$class_fwd[temp$aa_edited_fwd=="*"]     <- "nonsense"
#     temp$class_rev <- NA
#     temp$class_rev[temp$aa_edited_rev==temp$aa] <- "silent"
#     temp$class_rev[temp$aa_edited_rev!=temp$aa] <- "missense"
#     temp$class_rev[temp$aa_edited_rev=="*"]     <- "nonsense"
#     cols <- c("bp_rna_id", 
#         "aa_edited_fwd", "aa_edited_rev", 
#         "nuc_edited_fwd", "nuc_edited_rev",
#         "class_fwd", "class_rev"
#     )
#     temp <- temp[,cols]
#     wh <- match(txAnn$bp_rna_id, temp$bp_rna_id)
#     temp <- temp[wh,]
#     temp$bp_rna_id <- NULL
#     txAnn <- cbind(txAnn, temp)
#     rownames(txAnn) <- NULL
#     return(txAnn)
# }


# #' @export
# prepareTilingData <- function(txAnn,
#                               guidesAnn,
#                               lfcs, 
#                               tech=c("crisprko","crisprbe"),
#                               nuclease=c("Cas9", "Cas12a"),
#                               cbe_offset=-15,
#                               cbe_smooth=TRUE,
#                               cbe_type=c("c2t"),
#                               fun_aggregate=c("mean", "median"),
#                               remove_invalid=TRUE,
#                               aggregate_by_amino=TRUE,
#                               ...
# ){
#     tech     <- match.arg(tech)
#     nuclease <- match.arg(nuclease)
#     cbe_type    <- match.arg(cbe_type)
#     fun_aggregate <- match.arg(fun_aggregate)
#     txAnn  <- txAnn[!is.na(txAnn$aa_id),,drop=FALSE]
#     guides <- intersect(guidesAnn$ID, names(lfcs))
#     guidesAnn <- guidesAnn[match(guides, guidesAnn$ID),,drop=FALSE]

#     if (tech=="crisprko"){
#         guidesAnn$action_site <- .getCutSite(pam_site=guidesAnn$pam_site,
#                                              strand=guidesAnn$strand, 
#                                              nuclease=nuclease)
#     } else if (tech=="crisprbe"){
#         guidesAnn$action_site <- .getCutSite(pam_site=guidesAnn$pam_site,
#                                              strand=guidesAnn$strand, 
#                                              nuclease=nuclease,
#                                              cut_offset=cbe_offset)
#     }   
#     if (tech=="crisprko" | (tech=="crisprbe" & !cbe_smooth)){
#         guidesAnn <- guidesAnn[guidesAnn$action_site %in% txAnn$pos,,drop=FALSE]
#         df <- txAnn
#         wh <- match(df$pos, guidesAnn$action_site)
#         df$sgrna <- guidesAnn$ID[wh]
#         df$lfc <- lfcs[df$sgrna]
#     #Smoothing using kernel weights:
#     } else {
#         ws <- getCbeWeights()
#         nw <- (length(ws)-1)/2
#         smoothLfcs <- function(guides){
#             guides <- guides[!sapply(guides, is.null)]
#             guides <- lapply(guides, function(x){
#                 wh <- match(x$ID,names(lfcs))
#                 x$lfc <- lfcs[wh]
#                 x
#             })
#             lfcs   <- sapply(guides, function(x){
#                 if (nrow(x)==0){
#                     x <- NA
#                 } else {
#                     x <- sum(ws[as.character(x$dist)]*x$lfc, na.rm=TRUE)/sum(ws[as.character(x$dist)], na.rm=TRUE)
#                 }
#                 return(x)
#             })
#         }

#         guides <- lapply(txAnn$pos, function(x){
#             diff <- guidesAnn$action_site-x
#             wh <- which(abs(diff)<=nw)
#             if (length(wh)>0){
#                 temp <- data.frame(guidesAnn$ID[wh], diff[wh])
#                 colnames(temp) <- c("ID", "dist")
#             } else {
#                 temp <- NULL
#             }
#             return(temp)
#         })
#         names(guides) <- txAnn$bp_plot_id
#         guides <- guides[!sapply(guides, is.null)]
#         guides_fwd <- lapply(guides, function(x){
#             x$strand <- guidesAnn$strand[match(x$ID, guidesAnn$ID)]
#             x <- x[x$strand=="+",,drop=FALSE]
#             x
#         })
#         guides_rev <- lapply(guides, function(x){
#             x$strand <- guidesAnn$strand[match(x$ID, guidesAnn$ID)]
#             x <- x[x$strand=="-",,drop=FALSE]
#             x
#         })
#         lfcs_fwd <- smoothLfcs(guides_fwd)
#         lfcs_rev <- smoothLfcs(guides_rev)
#         wh1 <- match(names(lfcs_fwd), txAnn$bp_plot_id)
#         wh2 <- match(names(lfcs_rev), txAnn$bp_plot_id)
#         txAnn$lfc_fwd <- NA
#         txAnn$lfc_rev <- NA
#         txAnn$lfc_fwd[wh1] <- lfcs_fwd
#         txAnn$lfc_rev[wh2] <- lfcs_rev

#         if (cbe_type=="c2t" & remove_invalid){
#             txAnn$lfc_fwd[txAnn$nuc!="C"] <- NA
#             txAnn$lfc_rev[txAnn$nuc!="G"] <- NA
#         }
#         if (aggregate_by_amino){
#             dfs <- split(txAnn[,c("lfc_fwd","lfc_rev")], f=txAnn$aa_id)
#             dfs <- sapply(dfs, function(x){
#                 x <- unlist(x)
#                 x <- x[is.finite(x)]
#                 if (fun_aggregate=="mean"){
#                     x <- mean(x, na.rm=TRUE)
#                 }
#                 if (fun_aggregate=="median"){
#                     x <- median(x, na.rm=TRUE)
#                 }
#                 x
#             })
#             dfs <- dfs[!is.na(dfs)]
#             df <- txAnn
#             df$lfc <- NA
#             df$lfc[match(names(dfs), df$aa_id)] <- dfs
#         } else {
#             df <- txAnn
#             df$lfc <- sapply(1:nrow(df), function(i){
#                 mean(unlist(df[i, c("lfc_fwd", "lfc_rev")]), na.rm=TRUE)
#             })
#             df$lfc[!is.finite(df$lfc)] <- NA
#         }
#     }
#     df <- df[!is.na(df$lfc),,drop=FALSE]
#     df
# }

# #' @export
# #' @importFrom stats loess
# #' @importFrom stats predict
# plotTilingData <- function(df, loess=TRUE, ...){
#     x  <- df$aa_id
#     y  <- df$lfc
#     plot(x,y, pch=20, cex=0.5, col="grey75",...)
#     if (loess){
#         model  <- loess(y~x, span=300/max(x, na.rm=TRUE))
#         fitted <- predict(model, x)
#         lines(x,fitted,col="red", lwd=3)
#     }
# }




# #' @title Get binary segmentation results
# #' @description Get binary segmentation results from y (logFC values) against x (position vector)
# #' 
# #' @param x Position values (integer vector)
# #' @param y LogFC values (numeric vector)
# #' @param smooth.region Smoothing parameter for CNA algorithm.
# #' @param alpha Significance threshold for CNA algorithm.
# #' 
# #' @return A data.frame containing segments coordinates
# #' 
# #' @importFrom DNAcopy CNA smooth.CNA segment
# #' @export
# getBinarySegmentation <- function(x,
#                                   y,
#                                   smooth.region=1,
#                                   alpha=0.01
# ){
    
#     good <- !is.na(x) & !is.na(y)
#     x <- x[good]
#     y <- y[good]
#     lfc <- y
#     A <- median(lfc, na.rm=TRUE)
#     lfc <- lfc-A
#     chr <- rep("chr1", length(lfc))
#     cna <- suppressWarnings(CNA(cbind(lfc),
#         chrom=chr,
#         maploc=x,
#         data.type="logratio")
#     )
#     cna <- smooth.CNA(cna, smooth.region=smooth.region)
#     seg <- segment(cna, verbose=0, alpha=alpha)
#     out <- seg$output
#     out$seg.mean <- out$seg.mean+A
#     out$ID <- out$chrom <- NULL
#     out$num.mark <- NULL
#     return(out)
# }


# #' @importFrom graphics segments
# #' @export
# drawSegs <- function(segs, lwd=3, lty=1, col="firebrick3", ...){
#     for (i in 1:nrow(segs)){
#         segments(x0=segs[i,"loc.start"],
#                  x1=segs[i,"loc.end"],
#                  y0=segs[i,"seg.mean"],
#                  col=col,
#                  lwd=lwd,
#                  lty=lty,
#                  ...
#         )
#     }
# }



# .getCutSite <- function(pam_site,
#                         strand,
#                         nuclease=c("Cas9", "Cas12a"),
#                         cut_offset=NULL
# ){
#     nuclease <- match.arg(nuclease)
#     stopifnot(length(pam_site)==length(strand))
#     if (is.null(cut_offset)){
#         cut_offset <- .getDefaultCutOffset(nuclease)
#     } 
#     cut_site <- pam_site
#     cut_site[strand=='+'] <- pam_site[strand=='+'] + cut_offset
#     cut_site[strand=='-'] <- pam_site[strand=='-'] - cut_offset
#     return(cut_site)
# }




# .getDefaultCutOffset <- function(nuclease=c("Cas9", "Cas12a")){
#     nuclease <- match.arg(nuclease)
#     if (nuclease=="Cas9"){
#         offset=-4
#     } else if (nuclease=="Cas12a"){
#         offset=22
#     }
#     return(offset)
# }







#tab <- getTxTable("ENST00000000233")
#plot(tab$bp_plot_id, rep(1, nrow(tab)), pch=20, cex=0.1)
#wh <- which(!is.na(tab$bp_rna_id))
#abline(v=wh)
#wh <- which(!is.na(tab$aa))
#abline(v=wh, col="red")



