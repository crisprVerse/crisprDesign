#' @title Add gene context annotation to a \linkS4class{GuideSet} object
#' @description Add gene context annotation to spacer sequence stored
#'     in a \linkS4class{GuideSet} object
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param txObject A \linkS4class{TxDb} object or a
#'     \linkS4class{GRangesList} object obtained using
#'     \code{\link{TxDb2GRangesList}} to provide a 
#'     gene model annotation. 
#' @param anchor String specifying which relative coordinate
#'     of gRNAs should be used to locate gRNAs within gene.
#'     Must be either "cut_site", "pam_site" or "editing_site".
#' @param ignore_introns Should gene introns be ignored when annotating?
#'     TRUE by default. 
#' @param ignore.strand Should gene strand be ignored when annotating?
#'     TRUE by default. 
#' @param addPfam Should Pfam domains annotation be added?
#'     FALSE by default. If set to TRUE, \pkg{biomaRt} must be installed.
#' @param mart_dataset String specifying dataset to be used by \pkg{biomaRt}
#'     for Pfam domains annotation . E.g. "hsapiens_gene_ensembl".
#' @param ... Additional arguments, currently ignored.
#' 
#' @return A \linkS4class{GuideSet} object with a "geneAnnotation" list column
#'     stored in \code{mcols(guideSet)}. See details section for a
#'     description of the different gene annotation columns.
#' 
#' 
#' @details 
#' 
#' For DNA-targeting nucleases, the different columns stored in
#' \code{mcols(guideSet)[["geneAnnotation"]]} are:
#' 
#' \itemize{
#' \item \code{tx_id} Transcript ID.
#' \item \code{gene_symbol} Gene symbol.
#' \item \code{gene_id} Gene ID.
#' \item \code{protein_id} Protein ID.
#' \item \code{ID} gRNA ID.
#' \item \code{pam_site} gRNA PAM site coordinate.
#' \item \code{cut_site} gRNA cut site coordinate.
#' \item \code{chr} gRNA chromosome name.
#' \item \code{strand} gRNA strand. 
#' \item \code{cut_cds} Is the gRNA cut site located within the coding sequence
#'     (CDS) of the targeted isoform?
#' \item \code{cut_fiveUTRs} Is the gRNA cut site located within the 5'UTR
#'     of the targeted isoform?
#' \item \code{cut_threeUTRs} Is the gRNA cut site located within the 3'UTR
#'     of the targeted isoform?
#' \item \code{cut_introns} Is the gRNA cut site located within an intron
#'     of the targeted isoform?
#' \item \code{percentCDS} Numeric value to indicate the relative position of 
#'    the cut site with respect to the start of the CDS sequence when
#'    \code{cut_cds} is \code{TRUE}. The relative position is expressed as 
#'    a percentage from the total length of the CDS.
#' \item \code{percentTx} Numeric value to indicate the relative position of 
#'    the cut site with respect to the start of the mRNA sequence
#'    (therefore including 5' UTR). The relative position is expressed as a
#'    percentage from the total length of the mRNA sequence.
#' \item \code{aminoAcidIndex} If \code{cut_cds} is \code{TRUE}, integer value
#'     indicating the amino acid index with respect to the start of the protein.
#' \item \code{downstreamATG} Number of potential reinitiation sites 
#'     (ATG codons) downstream of the gRNA cut site, within 85 amino acids.
#' \item \code{nIsoforms} Numeric value indicating the number of isoforms
#'     targeted by the gRNA.
#' \item \code{totalIsoforms} Numeric value indicating the total number of
#'     isoforms existing for the gene targeted by the gRNA and specified
#'     in \code{gene_id}.
#' \item \code{percentIsoforms} Numeric value indicating the percentage of 
#'     isoforms for the gene specified in \code{gene_id} targeted by the gRNA.
#'     Equivalent to \code{nIsoforms}/\code{totalIsoforms}*100.
#' \item \code{isCommonExon} Logical value to indicate whether or not the gRNA
#'     is targeing an exon common to all isoforms.
#' \item \code{nCodingIsoforms} Numeric value indicating the number of 
#'     coding isoforms targeted by the gRNA. 5' UTRs and 3' UTRs are excluded.
#' \item \code{totalCodingIsoforms} Numeric value indicating the total number of
#'     coding isoforms existing for the gene targeted by the gRNA and specified
#'     in \code{gene_id}.
#' \item \code{percentCodingIsoforms} Numeric value indicating the percentage
#'     of coding isoforms for the gene specified in \code{gene_id} targeted
#'     by the gRNA. Equivalent to
#'     \code{nCodingIsoforms}/\code{totalCodingIsoforms}*100.
#'     5' UTRs and 3' UTRs are excluded.
#' \item \code{isCommonCodingExon} Logical value to indicate whether or
#'     not the gRNA is targeing an exon common to all coding isoforms.
#' }
#' 
#' @examples 
#' data(guideSetExample, package="crisprDesign")
#' data(grListExample, package="crisprDesign")
#' guideSet <- addGeneAnnotation(guideSetExample[1:6],
#'                               txObject=grListExample)
#' 
#' # To access a gene annotation already added:
#' ann <- geneAnnotation(guideSet)
#' 
#' @author Jean-Philippe Fortin, Luke Hoberecht
#' 
#' @seealso \code{\link{addTssAnnotation}} to add TSS annotation, and
#'     \code{\link{geneAnnotation}} to retrieve an existing gene annotation.
#' 
#' @export
#' @importFrom S4Vectors split mcols<-
#' @importFrom BiocGenerics rownames
#' @rdname addGeneAnnotation
setMethod("addGeneAnnotation",
          "GuideSet", 
          function(object,
                   txObject,
                   anchor=c("cut_site", "pam_site", "editing_site"),
                   ignore_introns=TRUE,
                   ignore.strand=TRUE,
                   addPfam=FALSE,
                   mart_dataset=NULL
){
    object <- .validateGuideSet(object)
    anchor <- match.arg(anchor)
    geneAnn <- .getGeneAnnotation(guideSet=object,
                                  txObject=txObject,
                                  anchor=anchor,
                                  ignore_introns=ignore_introns,
                                  ignore.strand=ignore.strand,
                                  addPfam=addPfam,
                                  mart_dataset=mart_dataset)
    splitFactor <- factor(BiocGenerics::rownames(geneAnn),
                          levels=names(object))
    geneAnn <- S4Vectors::split(geneAnn, f=splitFactor)
    S4Vectors::mcols(object)[["geneAnnotation"]] <- geneAnn
    return(object)
})



#' @rdname addGeneAnnotation
#' @export
setMethod("addGeneAnnotation",
          "PairedGuideSet", 
          function(object,
                   txObject,
                   anchor=c("cut_site", "pam_site", "editing_site"),
                   ignore_introns=TRUE,
                   ignore.strand=TRUE,
                   addPfam=FALSE,
                   mart_dataset=NULL
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addGeneAnnotation(unifiedGuideSet,
                                         txObject=txObject,
                                         anchor=anchor,
                                         ignore_introns=ignore_introns,
                                         ignore.strand=ignore.strand,
                                         addPfam=addPfam,
                                         mart_dataset=mart_dataset)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    
    return(out)
})


#' @rdname addSNPAnnotation
#' @export
setMethod("addGeneAnnotation", "NULL", function(object){
    return(NULL)
})



# Obtain a data.frame containing gene annotation from 
# a GuideSet object
#' @importFrom crisprBase isRnase
.getGeneAnnotation <- function(guideSet,
                               txObject,
                               anchor,
                               ignore_introns,
                               ignore.strand,
                               addPfam,
                               mart_dataset
){
    nuc <- crisprNuclease(guideSet)
    guideSet <- .dropNtcs(guideSet)
    if (crisprBase::isRnase(nuc)){
        geneAnn <- .getGeneAnnotation_rna_nuclease(guideSet=guideSet,
                                                   txObject=txObject)
    } else {
        geneAnn <- .getGeneAnnotation_dna_nuclease(guideSet=guideSet,
                                                   txObject=txObject,
                                                   anchor=anchor,
                                                   ignore_introns=ignore_introns,
                                                   ignore.strand=ignore.strand,
                                                   addPfam=addPfam,
                                                   mart_dataset=mart_dataset)
    }
    return(geneAnn)
} 




# To add gene annotation for RNA-targeting nucleases (e.g. CasRx)
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#' @importClassesFrom S4Vectors DataFrame
.getGeneAnnotation_rna_nuclease <- function(guideSet,
                                            txObject
){ 
    hasAlignments <- "alignments" %in% colnames(S4Vectors::mcols(guideSet))
    if (!hasAlignments){
        stop("For RNA-targeting nucleases, addSpacerAlignments has to be",
             " called before addGeneAnnotation.")
    }
    roster <- data.frame(tx=as.character(GenomeInfoDb::seqnames(guideSet)))
    key <- .getTx2GeneTable(txObject)
    if (any(!roster$tx %in% key$tx_id)){
        stop("Some transcripts stored in seqnames(guideSet) are not found in ",
             "the txObject.")
    }
    roster$gene_id <- key$gene_id[match(roster$tx, key$tx_id)]
    roster$ID <- names(guideSet)
    geneids <- unique(roster$gene_id)
    dfs <- split(roster, f=roster$gene_id)[geneids]
    aln <- onTargets(guideSet, unlist=TRUE)
    txTables <- lapply(dfs, function(df){
        geneid <- df$gene_id[1]
        txs <- unique(key$tx_id)[key$gene_id==geneid]
        guideSetSubset <- guideSet[df$ID]
        alnSubset <- aln[as.character(GenomeInfoDb::seqnames(aln)) %in% txs] 
        alnSubset <- cbind(ID=names(alnSubset),
                           tx=as.character(seqnames(alnSubset)))
        alnSubset <- as.data.frame(alnSubset)
        alns <- split(alnSubset, f=alnSubset$ID)
        ns <- vapply(alns, function(x){
            length(unique(x$tx))
        }, FUN.VALUE=numeric(1))
        txs <- vapply(alns, function(x) {
            paste0(unique(x$tx), collapse=",")
        }, FUN.VALUE=character(1))
        out <- data.frame(ID=names(txs),
                          targetedTxs=txs,
                          nTargetedTxs=ns)
        return(out)
        #ns <- table(data.frame(alnSubset)$ID)
        #return(ns)
    })
    names(txTables) <- NULL
    txTable <- do.call(rbind, txTables)
    txTable <- txTable[names(guideSet),]

    # Building the final data.frame
    txTable$tx_id <- roster$tx
    txTable$gene_id <- roster$gene_id
    
    dfs <- split(key, f=key$gene_id)
    txTable$nTotalTxs <- vapply(dfs, nrow, FUN.VALUE=numeric(1))[txTable$gene_id]
    txTable$percentTargetedTxs <- txTable$nTargetedTxs / txTable$nTotalTxs * 100
    
    # OK ready to add summary:
    cols <- c("ID", "tx_id", "gene_id",
              "targetedTxs", "nTargetedTxs", 
               "nTotalTxs", 
               "percentTargetedTxs")
    txTable <- txTable[,cols]
    txTable <- S4Vectors::DataFrame(txTable)
    return(txTable)
}




# To add gene annotation for DNA-targeting nucleases (SpCas9, enAsCas12a, etc)
#' @importFrom GenomeInfoDb seqlevelsStyle<- seqlevelsStyle
#' @importFrom S4Vectors isTRUEorFALSE
.getGeneAnnotation_dna_nuclease <- function(guideSet,
                                            txObject,
                                            anchor,
                                            ignore_introns,
                                            ignore.strand,
                                            addPfam,
                                            mart_dataset
){
    txObject <- .validateGRangesList(txObject)
    GenomeInfoDb::seqlevelsStyle(txObject) <-
        GenomeInfoDb::seqlevelsStyle(guideSet)
    stopifnot("'ignore.strand' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(ignore.strand)
    })

    if (targetOrigin(guideSet)=="customSequences"){
        stop("addGeneAnnotation is not available for custom sequences.")
    }
    bsgenome <- bsgenome(guideSet)
    
    geneAnn <- .annotateGeneOverlaps(guideSet=guideSet,
                                     txObject=txObject,
                                     anchor=anchor,
                                     ignore_introns=ignore_introns,
                                     ignore.strand=ignore.strand)
    geneAnn <- .addCutRegions(geneAnn=geneAnn,
                              txObject=txObject,
                              ignore.strand=ignore.strand)
    geneAnn <- .addCdsPositionAnnotation(geneAnn=geneAnn,
                                         txObject=txObject,
                                         bsgenome=bsgenome,
                                         ignore.strand=ignore.strand)
    geneAnn <- .addTxPositionAnnotation(geneAnn=geneAnn,
                                        bsgenome=bsgenome,
                                        txObject=txObject,
                                        ignore.strand=ignore.strand)
    geneAnn <- .addTranscriptIsoformSummary(geneAnn=geneAnn,
                                            txObject=txObject)
    geneAnn <- .addCodingIsoformSummary(geneAnn=geneAnn,
                                        txObject=txObject)
    geneAnn <- .addPfamDomains(geneAnn=geneAnn,
                               txObject=txObject,
                               addPfam=addPfam,
                               mart_dataset=mart_dataset)
    geneAnn <- .asDataFrame(geneAnn)
    return(geneAnn)
}



# Add annotation re. whether or not 
# the gRNAs cuts overlap a known gene
#' @importClassesFrom GenomicRanges GPos
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors mcols mcols<- isTRUEorFALSE queryHits subjectHits
.annotateGeneOverlaps <- function(guideSet,
                                  txObject,
                                  anchor,
                                  ignore_introns,
                                  ignore.strand
){
    anchor <- .validateAnchor(anchor, guideSet)
    anchorSites <- GenomicRanges::GPos(
        seqnames=GenomeInfoDb::seqnames(guideSet),
        pos=S4Vectors::mcols(guideSet)[[anchor]],
        strand=BiocGenerics::strand(guideSet))
    names(anchorSites)  <- names(guideSet)
    stopifnot("'ignore_introns' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(ignore_introns)
    })
    if (ignore_introns){
        targetAnnotation <- txObject[["exons"]]
    } else {
        targetAnnotation <- txObject[["transcripts"]]
    }
    overlaps <- GenomicRanges::findOverlaps(anchorSites,
                                            targetAnnotation,
                                            ignore.strand=ignore.strand)
    geneAnn <- anchorSites[S4Vectors::queryHits(overlaps)]
    targetAnnotationCols <- c("gene_symbol",
                              "gene_id",
                              "tx_id",
                              "protein_id",
                              "exon_id")
    for (i in targetAnnotationCols){
        indices <- S4Vectors::subjectHits(overlaps)
        S4Vectors::mcols(geneAnn)[[i]] <-
            S4Vectors::mcols(targetAnnotation)[[i]][indices]
    }
    return(geneAnn)
}




# Add annotation re. where in a gene the gRNA is cutting 
#' @importFrom S4Vectors mcols<-
.addCutRegions <- function(geneAnn,
                           txObject,
                           ignore.strand
){
    regions <- c('cds', 'fiveUTRs', 'threeUTRs', 'introns')
    for (i in regions){
        colname <- paste0('cut_', i)
        targetAnn <- txObject[[i]]
        S4Vectors::mcols(geneAnn)[[colname]] <- .spacersCutInRegion(
            geneAnn=geneAnn,
            targetAnn=targetAnn,
            ignore.strand=ignore.strand)
    }
    return(geneAnn)
}


# Helper function for .addCutRegions
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
.spacersCutInRegion <- function(geneAnn,
                                targetAnn,
                                ignore.strand
){
    overlaps <- GenomicRanges::findOverlaps(geneAnn,
                                            targetAnn,
                                            ignore.strand=ignore.strand)
    hasMatchingTxId <- geneAnn$tx_id[S4Vectors::queryHits(overlaps)] == 
        targetAnn$tx_id[S4Vectors::subjectHits(overlaps)]
    spacersCutInRegion <- vapply(seq_along(geneAnn), function(x){
        hits <- S4Vectors::queryHits(overlaps) == x
        any(hasMatchingTxId[hits])
    }, FUN.VALUE=logical(1))
    return(spacersCutInRegion)
}

# Add relative position within CDS where a gRNA cuts
#' @importFrom S4Vectors mcols<-
#' @importFrom GenomeInfoDb seqlengths
.addCdsPositionAnnotation <- function(geneAnn,
                                      txObject,
                                      bsgenome,
                                      ignore.strand
){
    cdsAnn <- .getTxAnnotationList(geneAnn=geneAnn,
                                   txObject=txObject,
                                   featureType="cds")
    aaSeq <- .getAminoAcidSequences(cdsAnn=cdsAnn,
                                    bsgenome=bsgenome)
    cdsPosAnn <- .getCdsPositionAnnotation(geneAnn=geneAnn,
                                           cdsAnn=cdsAnn,
                                           seqlengths=GenomeInfoDb::seqlengths(txObject),
                                           aaSeq=aaSeq,
                                           ignore.strand=ignore.strand)
    for (i in seq_along(cdsPosAnn)){
        mcolname <- names(cdsPosAnn)[i]
        S4Vectors::mcols(geneAnn)[[mcolname]] <- cdsPosAnn[[i]]
    }
    return(geneAnn)
}


# Add relative position within full mRNA where a gRNA cuts
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom BiocGenerics strand width
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges pos pintersect
#' @importClassesFrom IRanges IRanges
.addTxPositionAnnotation <- function(geneAnn,
                                     bsgenome,
                                     txObject,
                                     ignore.strand
){
    txAnn <- .getTxAnnotationList(geneAnn=geneAnn,
                                  txObject=txObject,
                                  featureType="exons")
    seqlengths <- GenomeInfoDb::seqlengths(txObject)
    percentTx <- rep(NA, length(geneAnn))
    inExons <- which(!is.na(S4Vectors::mcols(geneAnn)$exon_id))
    exonicGeneAnn <- geneAnn[inExons]
    
    txIds <- S4Vectors::mcols(geneAnn)$tx_id[inExons]
    txAnn <- txAnn[txIds]

    ## create GRanges encompassing all bases upstream of anchor sites
    chrs <- as.character(GenomeInfoDb::seqnames(exonicGeneAnn)) # identical w/ txAnn
    strands <- vapply(txAnn, function(x){
        as.character(unique(BiocGenerics::strand(x)))
    }, FUN.VALUE=character(1))
    limits <- seqlengths[chrs]
    limits[strands == "+"] <- 1
    allUpstream <- GenomicRanges::GRanges(
        seqnames=GenomeInfoDb::seqnames(exonicGeneAnn),
        ranges=IRanges::IRanges(
            start=pmin(limits, GenomicRanges::pos(exonicGeneAnn)),
            end=pmax(limits, GenomicRanges::pos(exonicGeneAnn))
        ),
        strand=BiocGenerics::strand(exonicGeneAnn)
    )
    
    ## find upstream bases in exons of target genes
    exonOverlaps <- GenomicRanges::pintersect(txAnn,
                                              allUpstream,
                                              ignore.strand=ignore.strand)
    txPos <- vapply(exonOverlaps, function(x){
        sum(BiocGenerics::width(x))
    }, FUN.VALUE=numeric(1))
    
    ## calculate percentTx
    txLens <- vapply(txAnn, function(x){
        sum(BiocGenerics::width(x))
    }, FUN.VALUE=numeric(1))
    percentTx[inExons] <- round(txPos / txLens * 100, 1)
    
    S4Vectors::mcols(geneAnn)[["percentTx"]] <- percentTx
    return(geneAnn)
}


# Get transcript annotation 
#' @importFrom S4Vectors mcols split
.getTxAnnotationList <- function(geneAnn,
                                 txObject,
                                 featureType
){
    txIds <- unique(geneAnn$tx_id)
    txIds <- txIds[!is.na(txIds)]
    txAnn <- queryTxObject(txObject=txObject,
                           featureType=featureType,
                           queryColumn="tx_id",
                           queryValue=txIds)
    # need to handle cases where length(txAnn)==0
    txAnn <- txAnn[order(S4Vectors::mcols(txAnn)$exon_rank)]
    ## remove superfluous mcols?
    txAnn  <- S4Vectors::split(txAnn, txAnn$tx_id)
    return(txAnn)
}



# Get amino sequences information
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings translate
#' @importClassesFrom Biostrings DNAStringSet
.getAminoAcidSequences <- function(cdsAnn,
                                   bsgenome
){
    nucleotideSequences <- BSgenome::getSeq(bsgenome, names=cdsAnn)
    nucleotideSequences <- lapply(nucleotideSequences, unlist)
    nucleotideSequences <- Biostrings::DNAStringSet(nucleotideSequences)
    aminoAcidSequences <- Biostrings::translate(nucleotideSequences)
    return(aminoAcidSequences)
}


# Get amino sequences information
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics strand width
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges pos pintersect
#' @importClassesFrom IRanges IRanges
#' @importFrom Biostrings subseq letterFrequency
.getCdsPositionAnnotation <- function(geneAnn,
                                      cdsAnn,
                                      seqlengths,
                                      aaSeq,
                                      ignore.strand
){
    percentCDS <- aminoAcidIndex <- downstreamATG <- rep(NA, length(geneAnn))
    ## get subset of annotations overlapping CDS
    inCDS <- which(S4Vectors::mcols(geneAnn)$cut_cds)
    geneAnn <- geneAnn[inCDS]
    txIds <- geneAnn$tx_id
    cdsAnn <- cdsAnn[txIds]
    aaSeq <- aaSeq[txIds]
    
    ## create GRanges encompassing all bases upstream of anchor sites
    chrs <- as.character(GenomeInfoDb::seqnames(geneAnn)) # identical w/ cdsAnn
    strands <- vapply(cdsAnn, function(x){
        as.character(unique(BiocGenerics::strand(x)))
    }, FUN.VALUE=character(1))
    limits <- seqlengths[chrs]
    limits[strands == "+"] <- 1
    allUpstream <- GenomicRanges::GRanges(
        seqnames=GenomeInfoDb::seqnames(geneAnn),
        ranges=IRanges::IRanges(
            start=pmin(limits, GenomicRanges::pos(geneAnn)),
            end=pmax(limits, GenomicRanges::pos(geneAnn))
        ),
        strand=BiocGenerics::strand(geneAnn)
    )
    
    ## find upstream bases in CDS of target genes
    cdsOverlaps <- GenomicRanges::pintersect(cdsAnn,
                                             allUpstream,
                                             ignore.strand=ignore.strand)
    txPos <- vapply(cdsOverlaps, function(x){
        sum(BiocGenerics::width(x))
    }, FUN.VALUE=numeric(1))
    
    ## calculate percentCDS, aminoAcidIndex
    cdsLens <- vapply(cdsAnn, function(x){
        sum(BiocGenerics::width(x))
    }, FUN.VALUE=numeric(1))
    percentCDS[inCDS] <- round(txPos / cdsLens * 100, 1)
    aminoAcidIndex[inCDS] <- ceiling(txPos / 3)
    
    ## calculate downstreamATG
    starts <- pmin(aminoAcidIndex[inCDS]+1, floor(cdsLens/3))
    ends <- pmin(aminoAcidIndex[inCDS]+84, floor(cdsLens/3)) # 84aa ~= 250nt
    aaSeq <- Biostrings::subseq(aaSeq,
                                start=starts,
                                end=ends)
    nMet <- Biostrings::letterFrequency(aaSeq, letters="M")
    nMet <- as.numeric(nMet)
    downstreamATG[inCDS] <- nMet
    
    cdsPositionAnnotation <- list("percentCDS"=percentCDS,
                                  "aminoAcidIndex"=aminoAcidIndex,
                                  "downstreamATG"=downstreamATG)
    return(cdsPositionAnnotation)
}






#' @importFrom S4Vectors mcols<-
.addTranscriptIsoformSummary <- function(geneAnn,
                                         txObject
){
    isoformAnn <- .getIsoformAnnotation(geneAnn=geneAnn,
                                        txObject=txObject,
                                        featureType="exons")
    for (i in colnames(isoformAnn)){
        S4Vectors::mcols(geneAnn)[[i]] <- isoformAnn[[i]]
    }
    return(geneAnn)
}


# Get a summary of how many protein-coding isoforms are targeted by a gRNA
#' @importFrom S4Vectors mcols<-
.addCodingIsoformSummary <- function(geneAnn,
                                     txObject
){
    isoformAnn <- .getIsoformAnnotation(geneAnn=geneAnn,
                                        txObject=txObject,
                                        featureType="cds")
    colnames(isoformAnn) <- gsub("Isoforms$", "CodingIsoforms",
                                 colnames(isoformAnn))
    colnames(isoformAnn) <- gsub("Exon$", "CodingExon",
                                 colnames(isoformAnn))
    for (i in colnames(isoformAnn)){
        S4Vectors::mcols(geneAnn)[[i]] <- isoformAnn[[i]]
    }
    return(geneAnn)
}


# Helper function for .addCodingIsoformSummary
#' @importFrom S4Vectors mcols
.getIsoformAnnotation <- function(geneAnn,
                                  txObject,
                                  featureType
){
    geneIds <- unique(geneAnn$gene_id)
    geneIds <- geneIds[!is.na(geneIds)]
    regions <- queryTxObject(txObject=txObject,
                             featureType=featureType,
                             queryColumn="gene_id",
                             queryValue=geneIds)
    
    isoformCountByGene <- S4Vectors::mcols(regions)[, c("gene_id", "tx_id")]
    isoformCountByGene <- unique(isoformCountByGene)
    isoformCountByGene <- table(isoformCountByGene$gene_id)
    totalIsoforms <- as.numeric(isoformCountByGene[geneAnn$gene_id])
    
    txTable <- geneAnn
    if (featureType=="cds"){
        txTable <- txTable[txTable$cut_cds]
    } 
    ids <- factor(names(txTable),
                  levels=unique(names(geneAnn)))
    nIsoforms <- table(spacer_id=ids,
                       gene_id=S4Vectors::mcols(txTable)[["gene_id"]])
    nIsoforms <- as.data.frame(nIsoforms)
    geneAnnInteraction <- interaction(names(geneAnn),
                                      S4Vectors::mcols(geneAnn)[["gene_id"]])
    nIsoformsInteraction <- interaction(nIsoforms$spacer_id,
                                        nIsoforms$gene_id)
    indexMatch <- match(geneAnnInteraction, nIsoformsInteraction)
    nIsoforms <- nIsoforms$Freq[indexMatch]
    
    percentIsoforms <- round(nIsoforms/totalIsoforms*100, 1)
    isCommonExon <- nIsoforms == totalIsoforms
    isoformAnnotation <- data.frame(nIsoforms=nIsoforms,
                                    totalIsoforms=totalIsoforms,
                                    percentIsoforms=percentIsoforms,
                                    isCommonExon=isCommonExon)
    return(isoformAnnotation)
}


# Add Pfam domain annotation (whether or not a gRNA cuts in a known Pfam domain)
#' @importFrom S4Vectors isTRUEorFALSE mcols mcols<-
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges pos
.addPfamDomains <- function(geneAnn,
                            txObject,
                            addPfam,
                            mart_dataset
){
    stopifnot("'addPfam' must be TRUE or FALSE" = {
        S4Vectors::isTRUEorFALSE(addPfam)
    })
    if (!addPfam){
        return(geneAnn)
    }
    if (!requireNamespace("biomaRt", quietly=TRUE)){
        message("Please install the biomaRt package to add Pfam annotation.")
    }
    
    mart <- biomaRt::useMart('ensembl')
    availableDatasets <- biomaRt::listDatasets(mart)$dataset
    if (is.null(mart_dataset) || !mart_dataset %in% availableDatasets){
        stop("mart_dataset '", mart_dataset, "' is not valid. ",
             "Check available datasets with listDatasets(useMart('ensembl')).")
    }
    mart <- biomaRt::useDataset(mart_dataset, mart=mart)
    
    # Get bm of Pfam domains from anchor coordinates:
    #chr <- GenomeInfoDb::seqnames(geneAnn)
    #chr <- gsub('[^0-9]', '', as.character(chr))
    attributes <- c('ensembl_transcript_id', 'pfam', 'pfam_start', 'pfam_end')
    #filters    <- c('chromosome_name', 'start', 'end')
    #values     <- list(chr,
    #                   IRanges::pos(geneAnn),
    #                   IRanges::pos(geneAnn))
    values <- unique(geneAnn$tx_id)
    filters  <- c('ensembl_transcript_id')
    cat("[addGeneAnnotation] Obtaining Pfam domains from Ensembl \n")
    bm <- biomaRt::getBM(attributes=attributes,
                         filters=filters,
                         values=values,
                         mart=mart)
    bm <- unique(bm)
    pfam <- lapply(seq_len(nrow(bm)), function(x){
        txId <- S4Vectors::mcols(geneAnn)$tx_id
        aaIndex <- S4Vectors::mcols(geneAnn)$aminoAcidIndex
        hasMatchingTx <- !is.na(txId) & txId == bm$ensembl_transcript_id[x]
        isNotBeforeDomain <- !is.na(aaIndex) & aaIndex >= bm$pfam_start[x]
        isNotAfterDomain <- !is.na(aaIndex) & aaIndex <= bm$pfam_end[x]
        inPfamDomain <- hasMatchingTx & isNotBeforeDomain & isNotAfterDomain
        spacerPfamDomain <- rep(NA, length(geneAnn))
        spacerPfamDomain[inPfamDomain] <- bm$pfam[x]
        spacerPfamDomain
    })
    pfam <- as.data.frame(pfam, row.names = NULL)
    pfam <- apply(pfam, 1, function(x){
        domains <- x[!is.na(x)]
        domains <- unique(domains)
        if (length(domains) > 0){
            paste0(domains, collapse=';')
        } else {
            NA
        }
    })
    S4Vectors::mcols(geneAnn)[["pfam"]] <- pfam
    return(geneAnn)
}
