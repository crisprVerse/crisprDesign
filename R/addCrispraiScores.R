#' @title Add CRISPRa/CRISPRi on-target scores to
#'     a \linkS4class{GuideSet} object.
#' @description Add CRISPRa/CRISPRi on-target scores to a
#'    \linkS4class{GuideSet} object. Only available for SpCas9, and for 
#'    hg38 genome. Requires \pkg{crisprScore} package to be installed.
#' 
#' @param object A \linkS4class{GuideSet} object or a 
#'     \linkS4class{PairedGuideSet} object.
#' @param gr A \linkS4class{GRanges} object derived from \code{queryTss} used
#'     to produce the \code{guideSet} object.
#' @param tssObject  A \linkS4class{GRanges} object containing TSS coordinates
#'     and annotation. The following columns must be present:
#'     "ID", promoter", "tx_id" and "gene_symbol".
#' @param geneCol String specifying which column of the \code{tssObject} should
#'     be used for a unique gene identified. "gene_id" by default.     
#' @param modality String specifying which modality is used.
#'     Must be either "CRISPRi" or "CRISPRa". 
#' @param chromatinFiles Named character vector of length 3 specifying
#'     BigWig files containing chromatin accessibility data. See
#'     crisprScore vignette for more information.
#' @param fastaFile String specifying fasta file of the hg38 genome.
#' @param ... Additional arguments, currently ignored.
#' 
#' @return \code{guideSet} with an added column for the CRISPRai score.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @seealso \code{\link{addOnTargetScores}} to add other on-target scores.
#' 
#' @export
#' @importFrom crisprScore getCrispraiScores
#' @rdname addCrispraiScores
setMethod("addCrispraiScores", "GuideSet",
          function(object,
                   gr,
                   tssObject,
                   geneCol="gene_id",
                   modality=c("CRISPRi", "CRISPRa"),
                   chromatinFiles=NULL,
                   fastaFile=NULL
){
    crisprNuclease <- crisprNuclease(object)
    data(SpCas9,
         package="crisprBase",
         envir=environment())
    if (!.identicalNucleases(crisprNuclease, SpCas9)){
        stop("[addCrispraiScores] Only SpCas9 is supported at the moment.")
    }

    if (genome(object)[1]!="hg38"){
        stop("addCrispraiScores] Only hg38 genome supported at the moment.")
    }
    modality <- match.arg(modality)
    tssFrame  <- .prepareTssFrame(tssObject,
                                  geneCol=geneCol)
    grnaFrame <- .prepareGrnaFrame(object, gr)
    scores <- crisprScore::getCrispraiScores(sgrna_df=grnaFrame,
                                             tss_df=tssFrame,
                                             chromatinFiles=chromatinFiles,
                                             fastaFile=fastaFile,
                                             modality=modality)
    scores <- scores[match(names(object), rownames(scores)),1]
    if (modality=="CRISPRa"){
        mcols(object)$score_crispra <- scores
    } else {
        mcols(object)$score_crispri <- scores
    }
    return(object)
})




#' @rdname addCrispraiScores
#' @export
setMethod("addCrispraiScores", "PairedGuideSet",
          function(object,
                   gr,
                   tssObject,
                   geneCol="gene_id",
                   modality=c("CRISPRi", "CRISPRa"),
                   chromatinFiles=NULL,
                   fastaFile=NULL
){
    object <- .validatePairedGuideSet(object)
    unifiedGuideSet <- .pairedGuideSet2GuideSet(object)
    unifiedGuideSet <- addCrispraiScores(unifiedGuideSet,
                                         gr=gr,
                                         tssObject=tssObject,
                                         geneCol=geneCol,
                                         modality=modality,
                                         chromatinFiles=chromatinFiles,
                                         fastaFile=fastaFile)
    out <- .addColumnsFromUnifiedGuideSet(object,
                                          unifiedGuideSet)
    return(out)
})



#' @rdname addCrispraiScores
#' @export
setMethod("addCrispraiScores", "NULL", function(object){
    return(NULL)
})





# Prepare a data.frame containing TSS information necessary
# to the CRISPRai algorithm
# The following columns are necessary:
# ID, gene_symbol, promoter, tx_id
.prepareTssFrame <- function(tssObject,
                             geneCol="gene_id"
){
    tssObject <- as.data.frame(tssObject)

    cols <- c("ID", "promoter","tx_id", geneCol)
    if (!all(cols %in% colnames(tssObject))){
        choices <- setdiff(cols, colnames(tssObject))
        stop("The following columns are missing in the tssObject: \n \t",
             paste0(choices, collapse=", "),".")
    }

    out <- data.frame(tss_id=tssObject$ID,
                      gene_symbol=tssObject[[geneCol]],
                      promoter=tssObject$promoter,
                      transcripts=tssObject$tx_id,
                      position=tssObject$start,
                      strand=tssObject$strand,
                      chr=tssObject$seqnames)

    # Check if there are any missing values:
    if (sum(is.na(out$gene_symbol))>0){
        stop("gene_symbol has some missing values.")
    }
    if (sum(out$gene_symbol=="")>0){
        stop("gene_symbol has some empty values.")
    }
    if (sum(is.na(out$promoter))>0){
        stop("promoter has some missing values.")
    }
    if (sum(out$promoter=="")>0){
        stop("promoter has some empty values.")
    }
    if (sum(is.na(out$position))>0){
        stop("start has some missing values.")
    }
    if (sum(is.na(out$strand))>0){
        stop("strand has some missing values.")
    }
    if (sum(is.na(out$seqnames))>0){
        stop("strand has some missing values.")
    }

    # Checking for final compatibility:
    good <- all(out$tss_id==paste0(out$gene_symbol, "_", out$promoter))
    if (!good){
        stop("The ID does not seem to be of the form geneCol_promoter.")
    }

    return(out)
}


# Prepare a data.frame containing gRNA information necessary
# to the CRISPRai algorithm from a guideSet object and a GRanges
# object that was used as input to create the GuideSet object.
.prepareGrnaFrame <-function(guideSet, gr){

    len <- spacerLength(guideSet)
    seqs <- spacers(guideSet, as.character=TRUE)
    if (len!=19){
        if (len==20){
            seqs <- substr(seqs, 2,20)
        } else{
            stop("spacer length must be of length 19 or 20.")
        }
    } 

    cols <- c("ID")
    colnames <- colnames(mcols(gr))
    if (!all(cols %in% colnames)){
        choices <- setdiff(cols, colnames)
        stop("The following columns are missing in the gr object: \n \t",
             paste0(choices, collapse=", "),".")
    }


    ids <- gr$ID[match(guideSet$region, names(gr))]
    if (sum(is.na(ids))>0){
        stop("Some of the guideSet regions cannot be found in the gr object.")
    }
    out <- data.frame(grna_id=names(guideSet),
                      tss_id=ids, 
                      pam_site=pamSites(guideSet),
                      strand=as.character(strand(guideSet)),
                      spacer_19mer=seqs)
    return(out)
}




