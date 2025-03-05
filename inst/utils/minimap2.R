minimap2 <- function(fastq,
                     outfile=NULL,
                     outdir=NULL,
                     index=NULL,
                     bin.minimap2=NULL,
                     bin.samtools=NULL,
                     bam=TRUE,
                     verbose=TRUE
){
    params <- "-ax map-ont --splice -g2000 -G200k -A1 -B2 -O2,32 -E1,0 -C0 -z200 -un --junc-bonus=0 --splice-flank=no"
    if (is.null(index)){
        index <- "minimap2_indices/GRCh38.mmi"
    }
    if (is.null(bin.minimap2)){
        bin.minimap2 <- system.file(package="crisprBinaries",
                                    "binaries/minimap2",
                                    mustWork=TRUE)
    }
    if (is.null(bin.samtools)){
        bin.samtools <- system.file(package="crisprBinaries",
                                    "binaries/samtools",
                                    mustWork=TRUE)
    }
    if (!is.null(outdir) & !dir.exists(outdir)){
        dir.create(outdir)
    }

    # Alignment with minimap:
    .minimap2_aligner <- function(fastq){
        if (is.null(outfile)){
            outfile <- gsub(".fastq|.fastq.gz", ".sam", basename(fastq))   
        }
        if (!is.null(outdir)){
            outfile <- file.path(outdir, outfile)
        }
        cmd <- paste0(bin.minimap2, " ",
            params, " ",
            index, " ",
            fastq, " > ",
            outfile
        )
        system(cmd)
    }

    .minimap2_aligner(fastq)
    if (bam){
        bamify(outfile,
               bin.samtools=bin.samtools)
    }
    if (verbose){
        message(paste0("Done with ", fastq))
    }
    NULL
}



bamify <- function(samfile,
                   bin.samtools=NULL
){
    if (is.null(bin.samtools)){
        bin.samtools <- system.file(package="crisprBinaries",
                                    "binaries/samtools",
                                    mustWork=TRUE)
    }
    bamfile <- gsub(".sam", ".bam", samfile)
    sorted.bamfile <- gsub(".bam", "sorted.bam", bamfile)
    .samToBam(samfile, bin.samtools=bin.samtools)
    .sortBam(bamfile, bin.samtools=bin.samtools)
    .indexSortedBam(sorted.bamfile, bin.samtools=bin.samtools)
}



# Convert SAM to BAM:
.samToBam <- function(samfile,
                      bin.samtools
){
    bamfile <- gsub(".sam", ".bam", samfile)
    cmd <- paste0(bin.samtools, 
                  " view -S -b ", samfile," > ",bamfile)
    system(cmd)
}

# Sort BAM file:
.sortBam <- function(bamfile,
                     bin.samtools
){
    output <- gsub(".bam", ".sorted.bam", bamfile)
    cmd <- paste0(bin.samtools,
                  " sort ",bamfile," -o ",output)
    system(cmd)
}

# Index sorted BAM file:
.indexSortedBam <- function(bamfile,
                            bin.samtools
){
    cmd <- paste0(bin.samtools,
                  " index ", bamfile)
    system(cmd)
}





#' @importFrom rtracklayer export 
#' @export
writeBam <- function(reads, outfile){
    if (!grepl(".bam", outfile)){
        outfile <- paste0(outfile, ".bam")
    }
    export(reads, format="bam", con=outfile)
}


#' @importFrom utils read.table
getReadLengthDistribution <- function(fastq){
    infile <- paste0(tempfile(), ".gz")

    # Copying to temporary location
    cat("[getReadLengthDistribution] copying \n")
    system(paste0("cp ",fastq, " ", infile))

    # Gunzipping
    cat("[getReadLengthDistribution] Unzipping \n")
    system(paste0("gunzip ",infile))
    infile <- gsub(".gz", "", infile)


    cat("[getReadLengthDistribution] Calculating distribution \n")
    program <- "awk"
    cmd <- "'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'"
    cmd <- paste0(program, " ", cmd, " ", infile, " > ", tempfile)  
    cat("[getReadLengthDistribution] Calculating distribution \n")
    system(cmd)


    # Getting summary table of read length frequencies:
    summary <- read.table(tempfile, header=FALSE)
    colnames(summary) <- c("length", "count")
    return(summary)
}



#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @export
getBamReads <- function(bamfile,
                        primaryOnly=TRUE,
                        addChrToName=TRUE){
    ### Get all reads data first:
    bamflag   <- scanBamFlag(isPaired=FALSE)
    bamparam  <- ScanBamParam(bamflag,
                              what=c("seq","flag", "mapq"))
    reads  <- readGAlignments(bamfile,
                              use.names=TRUE,
                              param=bamparam)
    # Only keeping primary alignments:
    if (primaryOnly){
        flag <- as.data.frame(reads)$flag
        toRemove <- c(2048, 2304, 2320, 2064, 256, 272)
        reads <- reads[!flag %in% toRemove]
    }

    if (addChrToName){
        if (!grepl("chr", seqlevels(reads)[1])){
            seqlevels(reads) <- paste0("chr", seqlevels(reads))
        }
    }
    return(reads)
}


#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @export
overlapReads <- function(reads, gr){
    hits <- findOverlaps(reads, gr)
    if (length(hits) != 0) {
        reads <- reads[unique(queryHits(hits))]
    } else {
        reads <- reads[0]
    }
    return(reads)
}

#' @importFrom BiocGenerics start end 
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
overlapReadsSpanningGuides <- function(reads, guides){
    gr <- GRanges(guides$chr[1],
                  IRanges(start=min(guides$pam_site, na.rm=TRUE),
                          end=max(guides$pam_site, na.rm=TRUE)))
    cond1 <- start(reads) <= start(gr)
    cond2 <- end(reads) >= end(gr)
    cond3 <- as.character(seqnames(reads))==as.character(seqnames(gr))[1]
    reads <- reads[cond1 & cond2 & cond3]
    return(reads)
}

#' @export
addReads <- function(crisprOntExperiment){
    sampleSheet <- crisprOntExperiment[["sampleSheet"]]
    reads <- lapply(sampleSheet$fileName, getBamReads)
    names(reads) <- sampleSheet$sampleName
    crisprOntExperiment[["reads"]] <- reads
    sampleSheet$nreads <- vapply(reads, length, FUN.VALUE=1)
    crisprOntExperiment[["sampleSheet"]] <- sampleSheet
    return(crisprOntExperiment)
}
#' @export
getReads <- function(crisprOntExperiment){
    return(crisprOntExperiment[["reads"]])
}

#' @importFrom S4Vectors List
#' @importFrom S4Vectors mcols
#' @export
addOnTargetReads <- function(crisprOntExperiment,
                             ont_margin=50
){
    crisprOntExperiment[["readsOn"]] <- List()
    df <- crisprOntExperiment[["sampleSheet"]]
    df$nreads.on  <- NA
    df$nreads.off <- NA

    for (k in seq_len(nrow(df))){
        gs <- crisprOntExperiment[["ontGuideSets"]][[df$ont[k]]]
        aln <- crisprDesign::offTargets(gs)
        gr_on <- GRanges(seqnames(gs),
                         IRanges(start=gs$cut_site-ont_margin,
                                 end=gs$cut_site+ont_margin))
        gr_off <- GRanges(seqnames(aln),
                      IRanges(start=aln$cut_site-ont_margin,
                              end=aln$cut_site+ont_margin))
        readsOn  <- overlapReads(crisprOntExperiment$reads[[k]], gr_on)
        readsOff  <- overlapReads(crisprOntExperiment$reads[[k]], gr_off)
        df$nreads.on[k]  <- length(readsOn)
        df$nreads.off[k] <- length(readsOff)
        crisprOntExperiment[["readsOn"]][[k]] <- readsOn
    }
    df$perc.on  <- df$nreads.on/df$nreads*100
    df$perc.off <- df$nreads.off/df$nreads*100
    df$ratio.on.off <- df$nreads.on/df$nreads.off
    crisprOntExperiment[["sampleSheet"]] <- df
    names(crisprOntExperiment[["readsOn"]]) <- names(crisprOntExperiment[["reads"]])
    return(crisprOntExperiment)
}







#' @importFrom crisprDesign getSpacerAlignments
#' @importFrom crisprDesignAux DataFrameToGuideSet
#' @importFrom utils data
buildGuideSet <- function(sequences,
                          bsgenome,
                          index,
                          crisprNuclease=NULL
){
    if (is.null(crisprNuclease)){
        data("SpCas9",
             package="crisprBase",
             envir=environment())
        crisprNuclease <- SpCas9
    }
    aln <- getSpacerAlignments(sequences,
                               bsgenome=bsgenome,
                               aligner_index=index,
                               n_mismatches=0)
    aln <- as.data.frame(aln)
    aln$ID <- rownames(aln)
    aln$width <- NULL
    aln$n_mismatches <- NULL
    aln$canonical <- NULL
    gs <- DataFrameToGuideSet(aln,
                              bsgenome=bsgenome,
                              crisprNuclease=crisprNuclease,
                              spacerCol="spacer")
    gs <- gs[order(start(gs))]
    names(gs) <- paste0("grna", 1:length(gs))
    return(gs)
}




#' @export
addKOEfficiencies <- function(crisprOntExperiment,
                              cut_margin=5
){
    df <- crisprOntExperiment[["sampleSheet"]]
    df$koEfficiencyProb <- NA
    for (k in seq_len(nrow(df))){
        if (!is.na(df[k,"controlSample"])){
            controlIndex <- which(df$sampleName==df[k, "controlSample"])
            reads_ko <- crisprOntExperiment[["readsOn"]][[k]]
            reads_wt <- crisprOntExperiment[["readsOn"]][[controlIndex]]
            guideSet <- crisprOntExperiment[["koGuideSets"]][[df[k,"ko"]]]
            ko <- calculateKOEfficiency(reads_wt=reads_wt,
                                        reads_ko=reads_ko,
                                        guideSet=guideSet,
                                        cut_margin=cut_margin)
            df$koEfficiencyProb[k] <- ko
        }
        print(k)
    }
    crisprOntExperiment[["sampleSheet"]] <- df
    return(crisprOntExperiment)
}


#' @importFrom GenomicAlignments sequenceLayer
#' @importFrom GenomicAlignments cigar qwidth
#' @importFrom S4Vectors DataFrame
#' @export
calculateKOEfficiency <- function(reads_wt,
                                  reads_ko,
                                  guideSet,
                                  keepInformativeOnly=TRUE,
                                  cut_margin
){
    ann <- DataFrame(id=names(reads_ko))
    ann$readLength <- qwidth(reads_ko)
    seqs_wt <- sapply(seq_along(reads_wt), function(i){
        as.character(sequenceLayer(mcols(reads_wt)$seq[i],
                                   cigar(reads_wt[i])))
    })
    seqs_ko <- sapply(seq_along(reads_ko), function(i){
        as.character(sequenceLayer(mcols(reads_ko)$seq[i],
                                   cigar(reads_ko[i])))
    })
    names(seqs_wt) <- names(reads_wt)
    names(seqs_ko) <- names(reads_ko)
    cut_sites <- guideSet$cut_site
    names(cut_sites) <- names(guideSet)

    levels <- c("A", "C","G", "T", ".", "-")
    .getWtProbabilities <- function(reads_wt,
                                    seqs_wt,
                                    cut_site){
        start <- cut_site - cut_margin
        end   <- cut_site + cut_margin
        temp <- data.frame(id=names(seqs_wt),
                           seq=seqs_wt,
                           start=start(reads_wt),
                           end=end(reads_wt))
        good <- vapply(seq_len(nrow(temp)), function(kk){
            temp[kk, "start"]<=start & temp[kk, "end"]>=end
        }, FUN.VALUE=TRUE)
        temp <- temp[good,,drop=FALSE]
        subseqs <- vapply(seq_len(nrow(temp)), function(kk){
            ss <- start - temp[kk, "start"]
            ee <- end - temp[kk, "start"]
            substr(temp[kk,"seq"],ss,ee)
        }, FUN.VALUE="a")
        probs <- list()
        len <- nchar(subseqs)[1]
        for (pos in seq_len(len)){
            x <- factor(substr(subseqs,pos,pos),levels=levels)
            x <- table(x)
            class(x) <- "numeric"
            probs[[pos]] <- x/sum(x)
        }
        probs <- do.call(cbind, probs)
        colnames(probs) <- seq_len(len)
        return(probs)
    }


    #cut_site <- cut_sites[[1]]
    .getCallPerCutSite <- function(cut_site,
                                   spacer_name){
        start <- cut_site - cut_margin
        end   <- cut_site + cut_margin
        pwt_control <- .getWtProbabilities(reads_wt, seqs_wt, cut_site)
        pko_control <- 1-pwt_control

        temp <- data.frame(id=names(seqs_ko),
                           seq=seqs_ko,
                           start=start(reads_ko),
                           end=end(reads_ko))
        good <- vapply(seq_len(nrow(temp)), function(kk){
            temp[kk, "start"]<=start & temp[kk, "end"]>=end
        }, FUN.VALUE=TRUE)
        temp <- temp[good,,drop=FALSE]
        subseqs <- vapply(seq_len(nrow(temp)), function(kk){
            ss <- start - temp[kk, "start"]
            ee <- end - temp[kk, "start"]
            substr(temp[kk,"seq"],ss,ee)
        }, FUN.VALUE="a")
        pwt <- list()
        pko <- list()
        len <- nchar(subseqs)[1]
        for (pos in seq_len(len)){
            x <- factor(substr(subseqs,pos,pos),levels=levels)
            pwt[[pos]] <- pwt_control[x,pos]
            pko[[pos]] <- pko_control[x,pos]
        }
        pwt <- do.call(cbind, pwt)
        pko <- do.call(cbind, pko)
        pwt <- apply(pwt,1,prod)
        pko <- apply(pko,1,prod)
        temp$KO <- pko>pwt
        temp <- temp[,c("id","KO")]
        wh <- match(ann$id, temp$id)
        colname <- paste0("ko_", spacer_name)
        ann[[colname]] <- temp$KO[wh]     
        return(ann)
    }
    for (cut_site_index in seq_along(cut_sites)){
        ann <- .getCallPerCutSite(cut_sites[[cut_site_index]],
                                  names(cut_sites)[[cut_site_index]])
    }
    ko_cols <- paste0("ko_", names(cut_sites))
    deletions <- as.matrix(ann[,ko_cols])
    informative <- rowSums(is.na(deletions))!=ncol(deletions)
    if (keepInformativeOnly){
        ann <- ann[informative,,drop=FALSE]
    }
    deletions <- as.matrix(ann[,ko_cols])
    ko <- rowSums(deletions, na.rm=TRUE)>=1
    out <- sum(ko)/length(ko)
    return(out)
}
















