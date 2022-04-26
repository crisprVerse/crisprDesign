#' @title To add edited alleles for a CRISPR base editing GuideSet
#' @description To add edited alleles for a CRISPR base editing GuideSet.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param baseEditor A \linkS4class{BaseEditor} object.
#' @param editingWindow A numeric vector of length 2 specifying
#'     start and end positions of the editing window with 
#'     respect to the PAM site. If \code{NULL} (default),
#'     the editing window of the \code{BaseEditor} object
#'     will be considered. 
#' @param nMaxAlleles Maximum number of edited alleles to report
#'     for each gRNA. Alleles from high to low scores.
#'     100 by default. 
#' @param addFunctionalConsequence Should variant classification
#'     of the edited alleles be added? TRUE by default.
#'     If \code{TRUE}, \code{txTable} must be provided.
#' @param addSummary Should a summary of the variant classified
#'     by added to the metadata columns of the \code{guideSet}
#'     object? TRUE by default. 
#' @param txTable Table of transcript-level nucleotide and amino
#'     acid information needed for variant classification.
#'     Usually returned by \code{\link{getTxInfoDataFrame}}.
#' @param verbose Should messages be printed to console?
#'     TRUE by default. 
#' 
#' @return The original \code{guideSet} object with an additional
#'     metadata column (\code{editedAlleles}) storing the annotated
#'     edited alelles. The edited alleles are always reported 
#'     from 5' to 3' direction on the strand corresponding to the
#'     gRNA strand. 
#' 
#' @examples
#' 
#' data(BE4max, package="crisprBase")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' txid="ENST00000357654"
#' library(crisprDesignData)
#' txTable <- getTxInfoDataFrame(tx_id=txid,
#'     txObject=txdb_human,
#'     bsgenome=bsgenome)
#' gr <- queryTxObject(txdb_human,
#'     queryValue=txid,
#'     queryColumn="tx_id",
#'     featureType="cds")
#' gs <- findSpacers(gr,
#'     bsgenome=bsgenome,
#'     crisprNuclease=BE4max)
#' gs <- unique(gs)
#' gs <- gs[1:2] # For the sake of time
#' gs <- addEditedAlleles(gs, baseEditor=BE4max,txTable=txTable)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
addEditedAlleles <- function(guideSet,
                             baseEditor,
                             editingWindow=NULL,
                             nMaxAlleles=100,
                             addFunctionalConsequence=TRUE,
                             addSummary=TRUE,
                             txTable=NULL,
                             verbose=TRUE
){
    if (addFunctionalConsequence & is.null(txTable)){
        stop("txTable must be provided when ",
             "addFunctionalConsequence is TRUE.")
    }

    if (verbose){
        message("[addEditedAlleles] Obtaining edited alleles at each gRNA target site.")
    }
    alleles <- lapply(seq_along(guideSet), function(guide){
        .getEditedAllelesPerGuide(gs=guideSet[guide],
                                  baseEditor=baseEditor,
                                  editingWindow=editingWindow,
                                  nMaxAlleles=nMaxAlleles)
    })
    if (addFunctionalConsequence){
        if (verbose){
            message("[addEditedAlleles] Adding functional consequences to alleles.")
        }
        alleles <- lapply(alleles,
                          .addFunctionalConsequences,
                          txTable)
    }
    mcols(guideSet)[["editedAlleles"]] <- alleles

    if (addSummary){
        guideSet <- .addSummaryFromEditingAlleles(guideSet)
    }
    return(guideSet)
}

.addSummaryFromEditingAlleles <- function(guideSet){
    alleles <- mcols(guideSet)[["editedAlleles"]]
    choices <- c("missense",
                 "nonsense",
                 "silent",
                 "not_targeting")
    scores <- lapply(alleles, function(x){
        out <- c(missense=0,
                 nonsense=0,
                 silent=0,
                 not_targeting=0)
        x <- split(x$score, f=x$variant)
        x <- vapply(x, sum, FUN.VALUE=1)
        out[names(x)] <- x
        return(out)
    })
    scores <- do.call(rbind, scores)
    scores <- scores[, seq_len(3), drop=FALSE]
    colnames(scores) <- paste0("score_", colnames(scores))
    variants <- .voteVariant(scores)
    mcols(guideSet)[colnames(scores)] <- scores
    mcols(guideSet)[["maxVariant"]] <- variants[["class"]]
    mcols(guideSet)[["maxVariantScore"]] <- variants[["score"]]
    return(guideSet)
}




.voteVariant <- function(scores){
    classes <- colnames(scores)
    classes <- gsub("score_", "", classes)
    pos <- apply(scores, 1, which.max)
    maxes <- apply(scores, 1, max)
    classes <- classes[pos]
    sums <- rowSums(as.matrix(scores), na.rm=TRUE)    
    classes[which(sums==0)] <- "not_targeting" 
    maxes[which(sums==0)] <- NA
    return(list(class=classes,
                score=maxes))
}








#' @importFrom crisprBase editingStrand
.getEditedAllelesPerGuide <- function(gs,
                                      baseEditor,
                                      editingWindow=c(-20,-8),
                                      nMaxAlleles=100
){
    if (!is(baseEditor, "BaseEditor")){
        stop("baseEditor must be a BaseEditor object.")
    }
    if (length(gs)!=1){
        stop("gs must be a GuideSet of length 1.")
    }
    if (editingStrand(baseEditor)!="original"){
        stop("Only base editors that edit the original ",
             "strand are supported at the moment. ")
    }
    ws <- .getEditingWeights(baseEditor,
                             editingWindow)
    nucChanges <- .getPossibleNucChanges(ws)


    # Getting gRNA information:
    pamSite <- pamSites(gs)
    strand <- as.character(strand(gs))
    chr <- as.character(seqnames(gs))

    seq <- .getExtendedSequences(gs, 
                                 start=editingWindow[1],
                                 end=editingWindow[2])
    nucs <- strsplit(seq, split="")[[1]]
    pos <- seq(editingWindow[1],
               editingWindow[2])
    names(nucs) <- pos


    # Getting scores for the edited nucleotides:
    nucsReduced <- nucs[nucs %in% names(nucChanges)]
    nucsReduced <- nucsReduced[names(nucsReduced) %in% colnames(ws)]
    ws <- ws[,colnames(ws) %in% names(nucsReduced),drop=FALSE]
    choices <- lapply(nucsReduced, function(x){
        nucChanges[[x]]
    })
    sequences <- expand.grid(choices)
    seqEdited <- apply(sequences, 1, paste0, collapse="")
    scores <- .scoreEditedAlleles(sequences,
                                  nucsReduced,
                                  ws)
  
    # Only keeping scores passing a threshold:
    reducedEditedAlleles <- data.frame(seq=seqEdited,
                                       score=scores)
    o <- order(-reducedEditedAlleles$score)
    reducedEditedAlleles <- reducedEditedAlleles[o,,drop=FALSE]
    sequences <- sequences[o,,drop=FALSE]


    nMaxAlleles <- min(nMaxAlleles, nrow(sequences))
    good <- seq_len(nMaxAlleles)
    reducedEditedAlleles <- reducedEditedAlleles[good,,drop=FALSE]
    sequences <- sequences[good,,drop=FALSE]
    sequences <- as.matrix(sequences)

    # Reconstructing full sequences:
    fullSequences <- c()
    for (i in seq_len(nrow(sequences))){
        temp <- nucs
        temp[colnames(sequences)] <- as.character(sequences[i,])
        temp <- lapply(temp, as.character)
        fullSequences[i] <- paste0(as.character(unlist(temp)),collapse="")
    }
    editedAlleles <- data.frame(seq=fullSequences,
                                score=reducedEditedAlleles$score)
    editedAlleles <- editedAlleles[order(-editedAlleles$score),,drop=FALSE]
    rownames(editedAlleles) <- NULL
    editedAlleles <- editedAlleles[editedAlleles$seq!=seq,,drop=FALSE]
    editedAlleles <- DataFrame(editedAlleles)

    # Adding metadata:
    metadata(editedAlleles)$wildtypeAllele <- seq
    if (strand=="+"){
        start <- pamSite + editingWindow[1]
        end   <- pamSite + editingWindow[2]
    } else {
        start <- pamSite - editingWindow[2]
        end   <- pamSite - editingWindow[1]
    }
    names(start) <- names(end) <- NULL
    metadata(editedAlleles)$start <- start
    metadata(editedAlleles)$end <- end
    metadata(editedAlleles)$chr <- chr
    metadata(editedAlleles)$strand <- strand
    metadata(editedAlleles)$editingWindow <- editingWindow
    editedAlleles$seq <- DNAStringSet(editedAlleles$seq)
    return(editedAlleles)
}




#' @importFrom Biostrings complement reverse
.addFunctionalConsequences <- function(editedAlleles,
                                       txTable
){
    if (txTable$chr[[1]]!=metadata(editedAlleles)$chr){
        stop("editedAlleles are not on the same chromosome.")
    }
    editedAlleles$variant <- "not_targeting"
    txTable <- txTable[txTable$region=="CDS",,drop=FALSE]
    geneStrand  <- metadata(txTable)$gene_strand
    guideStrand <- metadata(editedAlleles)$strand
    start <- metadata(editedAlleles)$start
    end <- metadata(editedAlleles)$end
    editingPositions <- start:end
    overlapPositions <- editingPositions[editingPositions %in% txTable$pos]

    if (length(overlapPositions)==0){
        return(editedAlleles)
    }

    # Getting nucleotide to replace
    sequences <- editedAlleles$seq
    if (geneStrand!=guideStrand){
        sequences <- complement(sequences)
    }
    if (guideStrand=="-"){
        sequences <- reverse(sequences)
    }
    nucs <- as.matrix(sequences)
    colnames(nucs) <- editingPositions
    nucs <- nucs[, as.character(overlapPositions), drop=FALSE]
    

    # Get wildtype protein:
    wh <- match(overlapPositions, txTable$pos)
    txTable <- txTable[order(txTable$pos_cds),,drop=FALSE]
    nuc <- txTable$nuc
    protein <- translate(DNAString(paste0(nuc, collapse="")))
    protein <- as.vector(protein)

    effects <- vapply(seq_len(nrow(nucs)), function(k){
        editedNuc <- nuc
        editedNuc[wh] <- nucs[k,]
        editedNuc <- DNAString(paste0(editedNuc, collapse=""))
        protein_edited <- as.vector(translate(editedNuc))

        mismatches <- which(protein_edited!=protein)
        if (length(mismatches)==0){
            effect <- "silent"
        } else {
            variants <- protein_edited[mismatches]
            if ("*" %in% variants){
                effect <- "nonsense"
            } else {
                effect <- "missense"
            }
        }
        return(effect)
    }, FUN.VALUE="a")

    aminos <- vapply(seq_len(nrow(nucs)), function(k){
        editedNuc <- nuc
        editedNuc[wh] <- nucs[k,]
        editedNuc <- DNAString(paste0(editedNuc, collapse=""))
        protein_edited <- as.vector(translate(editedNuc))
        # Getting amino acids
        aas <- rep(protein_edited, each=3)
        aas <- aas[wh]
        aas <- paste0(aas, collapse="")

        return(aas)
    }, FUN.VALUE="a")
    editedAlleles$variant <- effects
    editedAlleles$aa <- aminos

    # Adding wildtype amino:
    wildtypeAmino <- rep(protein, each=3)[wh]
    wildtypeAmino <- paste0(wildtypeAmino, collapse="")
    metadata(editedAlleles)$wildtypeAmino <- wildtypeAmino
    return(editedAlleles)
}












#' @importFrom crisprBase editingWeights
.getEditingWeights <- function(baseEditor,
                               editingWindow
){
    ws <- editingWeights(baseEditor)
    ws <- .rescaleWeights(ws)
    ws <- ws[, as.numeric(colnames(ws))>=editingWindow[1],drop=FALSE]
    ws <- ws[, as.numeric(colnames(ws))<=editingWindow[2],drop=FALSE]
    ws <- crisprBase:::.getReducedEditingMatrix(ws)
    ws <- .addWildtypeWeights(ws)
    return(ws)
}


.rescaleWeights <- function(ws){
    ws <- ws/max(ws, na.rm=TRUE)
    nucStart <- crisprBase:::.getOriginBaseFromRownames(rownames(ws))
    nucEnd   <- crisprBase:::.getTargetBaseFromRownames(rownames(ws))

    ind <- arrayInd(which.max(ws), c(nrow(ws),ncol(ws)))
    maxNuc <- nucStart[ind[1]]
    pos <- colnames(ws)[ind[2]]
    factor <- sum(ws[nucStart==maxNuc,pos])
    ws <- ws/factor
    return(ws)
}


.addWildtypeWeights <- function(ws){
    nucStart <- crisprBase:::.getOriginBaseFromRownames(rownames(ws))
    nucEnd   <- crisprBase:::.getTargetBaseFromRownames(rownames(ws))
    nucsStart <- unique(nucStart)
    addRows <- list()

    for (kk in seq_along(nucsStart)){
        nuc <- nucsStart[kk]
        wh <- which(nucStart==nuc)
        addRows[[kk]] <- 1-colSums(ws[wh,,drop=FALSE])
    }
    addRows <- do.call(rbind, addRows)
    rownames(addRows) <- paste0(nucsStart, "2", nucsStart)
    ws <- rbind(ws, addRows)
    return(ws)
}








.getPossibleNucChanges <- function(ws){
    ws_start <- crisprBase:::.getOriginBaseFromRownames(rownames(ws))
    ws_end   <- crisprBase:::.getTargetBaseFromRownames(rownames(ws))
    ws_pos <- as.integer(colnames(ws))
    nucs_that_can_changed <- unique(ws_start)
    nuc_choices <- lapply(nucs_that_can_changed, function(x){
        unique(c(ws_end[ws_start==x]), x)
    })
    names(nuc_choices) <- nucs_that_can_changed
    return(nuc_choices)
}





.scoreEditedAlleles <- function(sequences,
                                nucsReduced,
                                ws
){
    for (i in seq_len(ncol(sequences))){
        sequences[,i] <- paste0(nucsReduced[i],
                                "2",
                                sequences[,i])
    }
    scores <- matrix(0,
                     nrow=nrow(sequences),
                     ncol=ncol(sequences))
    for (i in seq_len(ncol(sequences))){
        pos <- colnames(sequences)[i]
        wh <- match(sequences[,i], rownames(ws))
        scores[,i] <- ws[wh,pos]
    }
    scores <- apply(scores,1, function(x){
        Reduce("*",x)
    })
    return(scores)
}









