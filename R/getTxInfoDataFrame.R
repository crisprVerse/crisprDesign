#' @title To obtain a DataFrame of transcript-specific CDS and mRNA coordinates
#' @description To obtain a DataFrame of transcript-specific CDS and mRNA
#'     coordinates.
#' 
#' @param tx_id String specifying ENSEMBL Transcript id.
#' @param txObject A \linkS4class{TxDb} object or a \linkS4class{GRangesList}
#'     object obtained using \code{\link{TxDb2GRangesList}}. 
#' @param bsgenome \linkS4class{BSgenome} object from which to extract
#'     sequences if a \linkS4class{GRanges} object is provided as input. 
#' @param extend Integer value specifying how many nucleotides in intron
#'     regions should be included.
#' @param checkCdsLength Should the CDS nucleotide length be a multiple of 3?
#'     TRUE by default. 
#' 
#' @return A \code{DataFrame} containing nucleotide and amino acid information.
#' The columns are:
#' 
#' \itemize{
#' \item \code{chr} Character specifying chromosome.
#' \item \code{pos} Integer value specifying coordinate in reference genome.
#' \item \code{strand} Character specifying strand of transcript.
#' \item \code{nuc} Character specifying nucleotide on the strand 
#'     specified by \code{strand}.
#' \item \code{aa} Character specifying amino acid.
#' \item \code{aa_number} Integer specifying amino acid number from 5' end.
#' \item \code{exon} Integer specifying exon number. 
#' \item \code{pos_plot} Integer specifying plot coordinate. Useful for
#'     plotting.
#' \item \code{pos_mrna} Integer specifying relative mRNA coordinate from the
#'     start of the mRNA.
#' \item \code{pos_cds} Integer specifying relative CDS coordinate from the
#'     start of the CDS.
#' \item \code{region} Character specifying gene region:
#'     3UTR, 5UTR, CDS, Intron, Upstream (promoter) or downstream. 
#' }
#' 
#' @examples 
#' 
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' data("grListExample")
#' tx_id <- "ENST00000538872"
#' df <- getTxInfoDataFrame(tx_id=tx_id,
#'     txObject=grListExample,
#'     bsgenome=bsgenome)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom Biostrings translate DNAString
#' @importFrom S4Vectors DataFrame
getTxInfoDataFrame <- function(tx_id,
                               txObject,
                               bsgenome,
                               extend=30,
                               checkCdsLength=TRUE
){  
    if (!tx_id %in% txObject$exons$tx_id){
        stop("transcript ID not found in txObject.")
    }
    if (extend<0){
        stop("extend must be a positive value.")
    }
    gr.exons <- txObject$exons[txObject$exons$tx_id==tx_id] 
    gr.exons <- gr.exons[order(gr.exons$exon_rank)]
    gr.cds   <- txObject$cds[txObject$cds$tx_id==tx_id]
    gr.cds   <- gr.cds[order(gr.cds$exon_rank)]
    if (!.isProteinCoding(gr.cds)){
        stop("The specified tx_id is not a protein coding transcript.")
    }
    if (checkCdsLength & !.hasValidProteinCodingLength(gr.cds)){
        stop("The specified tx_id has a CDS with incomplete length.")   
    }

    strand <- as.character(strand(gr.cds))[[1]]
    chr    <- as.character(seqnames(gr.cds))[[1]]


    df_aa <- .getCdsLongDataFrame(gr.cds=gr.cds,
                                  gr.exons=gr.exons,
                                  chr=chr,
                                  strand=strand,
                                  bsgenome=bsgenome)

    df_exon <- .getExonLongDataFrame(gr.exons=gr.exons,
                                     chr=chr,
                                     strand=strand,
                                     bsgenome=bsgenome)

    df <- merge(df_exon, df_aa, all=TRUE)
    if (extend>0){
        df_exon_extended <- .getExonLongDataFrameExtended(gr.exons=gr.exons,
                                                          chr=chr,
                                                          strand=strand,
                                                          bsgenome=bsgenome,
                                                          extend=extend)
        df <- merge(df, df_exon_extended, all=TRUE)
    } else {
        df[["pos_plot"]] <- df[["pos_mrna"]]
    }

    # Ordering by position:
    if (strand=="+"){
        df <- df[order(df$pos),,drop=FALSE]
    } else {
        df <- df[order(-df$pos),,drop=FALSE]
    }

    # Adding region:
    df$region <- "CDS"
    df$region[is.na(df$exon)] <- "Intron"
    df$region[df$region!="Intron" & is.na(df$aa)] <- "UTR"

    # Adding intergenic label:
    start_utr <- min(which(df$region=="UTR"))
    end_utr <- max(which(df$region=="UTR"))
    wh_intron <- which(df$region=="Intron")
    wh_upstream   <- wh_intron[wh_intron<start_utr]
    wh_downstream <- wh_intron[wh_intron>end_utr]
    df$region[wh_upstream]   <- "Upstream"
    df$region[wh_downstream] <- "Downstream"

    #Annotating UTR:
    start_cds <- min(which(df$region=="CDS"))
    df$region[df$region=="UTR" & seq_along(df$region)<start_cds] <- "5UTR"
    df$region[df$region=="UTR" & seq_along(df$region)>start_cds] <- "3UTR"

    # Reordering columns:
    cols <- c("chr", "pos", "strand",
              "nuc", "aa", "aa_number", "exon",  
              "pos_plot", "pos_mrna", "pos_cds",
              "region")
    df <- df[,cols,drop=FALSE]    
    
    # Reformatting
    out <- DataFrame(df)
    metadata(out)$tx_id <- tx_id
    metadata(out)$extend <- extend
    metadata(out)$bsgenome <- bsgenome@pkgname
    metadata(out)$gene_strand <- out$strand[1]
    out$strand <- NULL
    return(out)
}



# Make sure a given CDS is really coding by checking
# that the length is greater than 0
.isProteinCoding <- function(gr.cds){
    len <- sum(BiocGenerics::width(gr.cds))
    len>0
}

# Make sure the CDS has the proper length
# to produce a valid protein
.hasValidProteinCodingLength <- function(gr.cds){
    len <- sum(BiocGenerics::width(gr.cds))
    (len %% 3)==0   
}


# Get CDS annotation and genomic coordinates
.getCdsLongDataFrame <- function(gr.cds,
                                 gr.exons,
                                 chr,
                                 strand,
                                 bsgenome,
                                 posId="pos_cds"
){
    if (strand=="+"){
        pos <- lapply(seq_along(gr.cds), function(i){
            seq(start(gr.cds)[i], end(gr.cds)[i], 1)
        }) 
    } else {
        pos <- lapply(seq_along(gr.cds), function(i){
            rev(seq(start(gr.cds)[i], end(gr.cds)[i], 1))
        }) 
    }
    pos <- do.call("c", pos)
    out <- data.frame(chr=chr,
                      pos=pos,
                      strand=strand)
    out[[posId]] <- seq_along(pos)

    # Adding aa number:
    naa <- length(pos)/3
    aa  <- rep(seq_len(naa), each=3)
    out[["aa_number"]] <- aa
    out <- out[out$aa!=naa,,drop=FALSE]


    # Adding exon number:
    gr.temp <- GRanges(chr,
                       IRanges(start=out$pos,
                               end=out$pos))
    hits <- findOverlaps(gr.temp, gr.exons)
    if (length(hits)!=nrow(out)){
        stop("Exons and cds not compatible. ")
    }
    out$exon <- subjectHits(hits)

    # Getting CDS sequence:
    seq <- getSeq(bsgenome, gr.cds, as.character=TRUE)
    seq <- paste0(seq, collapse="")
    stopifnot(substr(seq,1,3)=="ATG")
    seq <- substr(seq, 1, nchar(seq)-3)
    
    # Translating:
    prot <- translate(DNAString(seq))
    prot <- as.character(prot)
    prot <- rep(strsplit(prot,"")[[1]], each=3)

    # Adding nucleotide to output:
    seq <- strsplit(seq, "")[[1]]
    out$nuc <- seq

    # Adding amino acid to output:
    out$aa <- prot
    cols  <- c("chr", "pos", "strand",
               "exon", "nuc", "aa", "aa_number", posId)
    out <- out[,cols,drop=FALSE]
    return(out)
}

# Get full mRNA (coding and not coding) annotation and genomic coordinates
.getExonLongDataFrame <- function(gr.exons,
                                  chr,
                                  strand,
                                  bsgenome,
                                  posId="pos_mrna"
){
    if (strand=="+"){
        pos <- lapply(seq_along(gr.exons), function(i) {
            seq(start(gr.exons)[i], end(gr.exons)[i], 1)
        }) 
    } else {
        pos <- lapply(seq_along(gr.exons), function(i) {
            rev(seq(start(gr.exons)[i], end(gr.exons)[i], 1))
        }) 
    }
    pos <- do.call("c", pos)
    pos <- pos[!duplicated(pos)] #When extend>0, removes overlapping introns:

    out <- data.frame(chr=chr,
                      pos=pos,
                      strand=strand)
    out[[posId]] <- seq_along(pos)

    # Adding exon number:
    gr.temp <- GRanges(chr,
                       IRanges(start=out$pos,
                               end=out$pos))
    hits <- findOverlaps(gr.temp, gr.exons)
    out$exon <- subjectHits(hits)
    
    # Adding sequence info:
    seq <- getSeq(bsgenome, gr.exons, as.character=TRUE)
    seq <- paste0(seq, collapse="")
    seq <- strsplit(seq, "")[[1]]
    out$nuc <- seq
    cols  <- c("chr", "pos", "strand", "nuc", "exon", posId)
    out <- out[, cols, drop=FALSE]
    return(out)
}

# Get full mRNA (coding and not coding) annotation and genomic coordinates
# with an extended windown
.getExonLongDataFrameExtended <- function(gr.exons,
                                          chr,
                                          strand,
                                          bsgenome,
                                          extend=30,
                                          posId="pos_plot"
){
    if (strand=="+") {
        pos <- lapply(seq_along(gr.exons), function(i){
            seq(start(gr.exons)[i]-extend,
                end(gr.exons)[i]+extend, 1)
        }) 
    } else {
        pos <- lapply(seq_along(gr.exons), function(i){
            rev(seq(start(gr.exons)[i]-extend,
                    end(gr.exons)[i]+extend, 1))
        }) 
    }
    pos <- do.call("c", pos)
    pos <- pos[!duplicated(pos)] # Removing overlaping introns
    out <- data.frame(chr=chr,
                      pos=pos,
                      strand=strand)
    out[[posId]] <- seq_along(pos)

    gr_extended <- GRanges(chr,
                           IRanges(start=out$pos,
                                   end=out$pos),
                           strand=out$strand)

    seq <- getSeq(bsgenome, gr_extended, as.character=TRUE)
    seq <- paste0(seq, collapse="")
    seq <- strsplit(seq, "")[[1]]
    out$nuc <- seq
    cols  <- c("chr", "pos", "strand", "nuc", posId)
    out <- out[, cols, drop=FALSE]
    return(out)
}
