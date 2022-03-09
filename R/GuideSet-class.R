#' An S4 class to store CRISPR gRNA sequences with modular annotations.
#' 
#' @section Constructors:
#'     Use the constructor \code{link{GuideSet}} to create a GuideSet object.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{crisprNuclease}:}{To get \linkS4class{CrisprNuclease} object
#'         used to design gRNAs.}
#'     \item{\code{spacers}:}{To get spacer sequences.}
#'     \item{\code{protospacers}:}{To get protospacer sequences.}
#'     \item{\code{spacerLength}:}{To get spacer length.}
#'     \item{\code{protospacerLength}:}{To get protospacer length.}
#'     \item{\code{pams}:}{To get PAM sequences.}
#'     \item{\code{pamSites}:}{To get PAM site coordinates.}
#'     \item{\code{pamLength}:}{To get PAM length.}
#'     \item{\code{pamSide}:}{To return the side of the PAM sequence with
#'         respect to the protospacer sequence.}
#'     \item{\code{prototypeSequence}:}{To get a prototype protospacer
#'         sequence.}
#'     \item{\code{cutSites}:}{To get cut sites.}
#'     \item{\code{alignments}:}{To get genomic alignments annotation.}
#'     \item{\code{onTargets}:}{To get on-target alignments annotation}
#'     \item{\code{offTargets}:}{To get off-target alignments annotation}
#'     \item{\code{snps}:}{Tp get SNP annotation.}
#'     \item{\code{geneAnnotation}:}{To get gene annotation.}
#'     \item{\code{tssAnnotation}:}{To get TSS annotation.}
#'     \item{\code{enzymeAnnotation}:}{To get restriction enzymes annotation.}
#' }
#' 
#' @export
setClass("GuideSet", contains = "GRanges")






#' @describeIn GuideSet Create a \linkS4class{GuideSet} object
#' @param protospacers Character vector of protospacers sequences.
#' @param pams Character vector of PAM sequences.
#' @param seqnames Character vector of chromosome names.
#' @param pam_site Integer vector of PAM site coordinates.
#' @param strand Character vector of gRNA strand.
#'    Only accepted values are "+" and "-".
#' @param CrisprNuclease \linkS4class{CrisprNuclease} object.
#' @param targetOrigin String specifying the origin of the DNA target.
#'     Must be either 'bsgenome' or 'customSequences'.
#' @param bsgenome \linkS4class{BSgenome} object. Must be specified when
#'     \code{targetOrigin} is equal to "bsgenome".
#' @param customSequences \linkS4class{DNAStringSet} object. Must be specified
#'     when \code{targetOrigin} is equal to "customSequences".
#' @param ... Additional arguments for class-specific methods
#' @param seqinfo A \linkS4class{Seqinfo} object containing informatioon
#'     about the set of genomic sequences present in the target genome.
#' @param seqlengths \code{NULL}, or an integer vector named with \code{levels(seqnames)}
#'     and containing the lengths (or NA) for each level in \code{levels(seqnames)}.
#' 
#' @return A GuideSet object.
#' @examples
#' protospacers <- c("AGGTCGTGTGTGGGGGGGGG",
#'                   "AGGTCGTGTGTGGGGGGGGG")
#' pams <- c("AGG", "CGG")
#' pam_site=c(10,11)
#' seqnames="chr7"
#' data(SpCas9, package="crisprBase")
#' CrisprNuclease <- SpCas9
#' strand=c("+", "-")
#' gr <- GuideSet(protospacers=protospacers,
#'                pams=pams,
#'                seqnames=seqnames,
#'                CrisprNuclease=CrisprNuclease,
#'                pam_site=pam_site,
#'                strand=strand,
#'                targetOrigin="customSequences",
#'                customSequences=protospacers)
#' 
#' @importFrom crisprBase CrisprNuclease
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom S4Vectors mcols<- mcols
#' @importFrom methods new
#' @export
GuideSet <- function(protospacers = NA_character_,
                     pams = NULL,
                     seqnames = NA_character_,
                     pam_site = 0L,
                     strand = "*",
                     CrisprNuclease = NULL,
                     targetOrigin = c("bsgenome", "customSequences"),
                     bsgenome = NULL,
                     customSequences = NULL,
                     ...,
                     seqinfo = NULL,
                     seqlengths = NULL
){
    targetOrigin <- match.arg(targetOrigin)
    protospacers <- .validateGuideSetSequences("protospacers", protospacers)
    pams <- .validateGuideSetSequences("pams", pams)
    gr <- GRanges(seqnames,
                  IRanges(start=pam_site,
                          width=1),
                  strand=strand,
                  ...,
                  seqinfo=seqinfo,
                  seqlengths=seqlengths)

    # Adding global metadata:
    metadata(gr)[["CrisprNuclease"]] <- CrisprNuclease
    metadata(gr)[["targetOrigin"]] <- targetOrigin
    if (targetOrigin=="bsgenome"){
        .isBSGenome(bsgenome)
        metadata(gr)[["bsgenome"]] <- bsgenome
    } else if (targetOrigin=="customSequences"){
        if (is.character(customSequences)){
            customSequences <- .string2DNAStringSet(customSequences,
                                                    both_strands=FALSE)
        }
        .isDNAStringSet(customSequences)
        metadata(gr)[["customSequences"]] <- customSequences
    }

    # Adding metadata columns:
    mcols(gr)[["protospacer"]] <- DNAStringSet(protospacers)
    if (!is.null(pams)){
        mcols(gr)$pam <- DNAStringSet(pams)
    }
    mcols(gr)[["pam_site"]] <- pam_site
    new("GuideSet", gr)
}



#' @importFrom crisprBase nucleaseName
#' @importFrom methods callNextMethod
setMethod("show",
          signature(object = "GuideSet"),
          function(object){
    callNextMethod()
    name <- nucleaseName(metadata(object)$CrisprNuclease)
    cat(paste0("  crisprNuclease: ", name, "\n"))
})



setValidity("GuideSet", function(object){

    df <- mcols(object)
    mandatoryCols <- c("protospacer", "pam_site","pam")
    out <- TRUE
    if (!all(mandatoryCols %in% colnames(df))){
        out <- paste0("The following columns must be present",
               " in mcols(object): protospacer, pam_site and pam.")
        return(out)
    } 

    if (!is(df[["protospacer"]], "DNAStringSet")){
        out <- "mcols(object)$protospacer must be a DNAStringSet object."
        return(out)
    } 
    if (!is(df[["pam"]], "DNAStringSet")){
        out <- "mcols(object)$pam must be a DNAStringSet object."
        return(out)
    } 

    meta <- metadata(object)
    mandatoryMetaFields <- c("CrisprNuclease", "targetOrigin")
    if (!all(mandatoryMetaFields %in% names(meta))){
        out <- paste0("The following field must be present",
                      " in metadata(object): CrisprNuclease and targetOrigin.")
        return(out)
    } 

    targetOriginChoices <- c("bsgenome", "customSequences")
    target <- meta[["targetOrigin"]]
    if (length(target)!=1){
        stop("targetOrigin must be a character vector of length 1.")
    }
    if (!target %in% targetOriginChoices){
        stop("targetOrigin must be either 'bsgenome' or 'customSequences'.")
    }

    if (!target%in% names(meta)){
        stop("When 'targetOrigin' is set to ", target, ", '",
         target, "' must be specified in the metadata field.")
    }    
    if (target=="bsgenome"){
        .isBSGenome(meta[["bsgenome"]])
    } else if (target=="customSequences"){
        .isDNAStringSet(meta[["customSequences"]])
    }


    nuc <- meta[["CrisprNuclease"]]
    if (!is(nuc, "CrisprNuclease")){
        out <- "metadata(object)$CrisprNuclease must be a CrisprNuclease object"
        return(out)
    }
    return(out)
})



#' @rdname GuideSet-class
#' @param object \linkS4class{GuideSet} object.
#' @export
setMethod("targetOrigin", "GuideSet", 
    function(object){
    out <- metadata(object)[["targetOrigin"]]
    return(out)
})



#' @rdname GuideSet-class
#' @param object \linkS4class{GuideSet} object.
#' @export
setMethod("customSequences", "GuideSet", 
    function(object){
    out <- metadata(object)[["customSequences"]]
    return(out)
})


#' @rdname GuideSet-class
#' @param object \linkS4class{GuideSet} object.
#' @export
setMethod("bsgenome", "GuideSet", 
    function(object){
    out <- metadata(object)[["bsgenome"]]
    return(out)
})












#' @rdname GuideSet-class
#' @param object \linkS4class{GuideSet} object.
#' @export
setMethod("crisprNuclease", "GuideSet", 
    function(object){
    out <- metadata(object)[["CrisprNuclease"]]
    return(out)
})


#' @rdname GuideSet-class
#' @param as.character Should sequences be returned as a character
#'     vector? FALSE by default, in which case sequences are returned
#'     as a \linkS4class{DNAStringSet}.
#' @param returnAsRna Should the sequences be returned as RNA
#'     instead of DNA? FALSE by default. 
#' @export
#' @importFrom crisprBase isRnase
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings RNAStringSet
setMethod("spacers", "GuideSet", 
    function(object,
             as.character=FALSE,
             returnAsRna=FALSE){
    out <- mcols(object)[["protospacer"]]
    if (returnAsRna){
        out <- RNAStringSet(out)
    }
    if (isRnase(crisprNuclease(object))){
        out <- reverseComplement(out)
    }
    if (as.character){
        out <- as.character(out)
    }
    names(out) <- names(object)
    return(out)
})



#' @rdname GuideSet-class
#' @export
setMethod("pams", "GuideSet", 
    function(object,
             as.character=FALSE,
             returnAsRna=FALSE){
    out <- mcols(object)[["pam"]]
    if (returnAsRna){
        out <- RNAStringSet(out)
    }
    if (as.character){
        out <- as.character(out)
    }
    names(out) <- names(object)
    return(out)
})

#' @rdname GuideSet-class
#' @export
setMethod("pamSites", "GuideSet", 
    function(object){
    out <- mcols(object)[["pam_site"]]
    names(out) <- names(object)
    return(out)
})



#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase cutSites
#' @importFrom BiocGenerics strand
setMethod("cutSites", "GuideSet", 
    function(object){
    pamSites <- mcols(object)[["pam_site"]]
    nuc <- metadata(object)[["CrisprNuclease"]]
    strand <- as.character(strand(object))
    out <- getCutSiteFromPamSite(pam_site=pamSites,
                                 strand=strand,
                                 crisprNuclease=nuc)
    names(out) <- names(object)
    return(out)
})





#' @rdname GuideSet-class
#' @param include.pam Should PAM sequences be included?
#'     FALSE by default. 
#' @export
setMethod("protospacers", "GuideSet", 
    function(object,
             as.character=FALSE,
             include.pam=FALSE,
             returnAsRna=FALSE){
    out <- mcols(object)[["protospacer"]]
    if (include.pam){
        pams <- pams(object)
        out  <- paste0(out, pams)
    }
    out <- DNAStringSet(out)
    if (returnAsRna){
        out <- RNAStringSet(out)
    }
    if (as.character){
        out <- as.character(out)
    }
    names(out) <- names(object)
    return(out)
})




#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase spacerLength
setMethod("spacerLength", "GuideSet", 
    function(object){
    nuc <- metadata(object)$CrisprNuclease
    out <- spacerLength(nuc)
    return(out)
})





#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase protospacerLength
setMethod("protospacerLength", "GuideSet", 
    function(object){
    nuc <- metadata(object)$CrisprNuclease
    out <- protospacerLength(nuc)
    return(out)
})

#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase prototypeSequence
setMethod("prototypeSequence", "GuideSet", 
    function(object){
    nuc <- metadata(object)$CrisprNuclease
    out <- prototypeSequence(nuc)
    return(out)
})



#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase pamLength
setMethod("pamLength", "GuideSet", 
    function(object){
    nuc <- metadata(object)$CrisprNuclease
    out <- pamLength(nuc)
    return(out)
})

#' @rdname GuideSet-class
#' @export
#' @importFrom crisprBase pamSide
setMethod("pamSide", "GuideSet", 
    function(object){
    nuc <- metadata(object)$CrisprNuclease
    out <- pamSide(nuc)
    return(out)
})






#' @rdname GuideSet-class
#' @param unlist Should the annotation be returned as
#'     one table instead of a list? TRUE by default. 
#' @importFrom BiocGenerics unlist
#' @export
setMethod("snps", "GuideSet", 
    function(object,
             unlist=TRUE){
    if (!"snps" %in% colnames(mcols(object))){
        out <- NULL
    } else {
        out <- mcols(object)[["snps"]]
        if (unlist){
            out <- BiocGenerics::unlist(out)
        }
    }
    return(out)
})


#' @rdname GuideSet-class
#' @param columnName Name of the column storing the alignments annotation
#'     to be retrieved.
#' @importFrom BiocGenerics unlist
#' @export
setMethod("alignments", "GuideSet", 
    function(object,
             columnName="alignments",
             unlist=TRUE){
    if (!columnName %in% colnames(mcols(object)) ||
        !.isAlignmentsColumn(object, columnName)){
        out <- NULL
    } else {
        out <- mcols(object)[[columnName]]
        if (unlist){
            out <- BiocGenerics::unlist(out)
            names(out) <- NULL
        }
    }
    return(out)
})



#' @rdname GuideSet-class
#' @importFrom BiocGenerics unlist
#' @export
setMethod("onTargets", "GuideSet", 
    function(object,
             columnName="alignments",
             unlist=TRUE){
    if (!columnName %in% colnames(mcols(object)) ||
        !.isAlignmentsColumn(object, columnName)){
        out <- NULL
    } else {
        out <- mcols(object)[[columnName]]
        out <- BiocGenerics::unlist(out)
        out <- out[out$n_mismatches == 0]
        names(out) <- NULL
        if (!unlist){
            spacers <- spacers(object, as.character=TRUE)
            out <- split(out, f=factor(out$spacer,
                                       levels=unique(spacers)))
            out <- out[spacers]
            names(out) <- names(object)
        }
    }
    return(out)
})


#' @rdname GuideSet-class
#' @param max_mismatches What should be the maximum number of 
#'     mismatches considered for off-targets? 
#'     Inf by default.
#' @importFrom BiocGenerics unlist
#' @export
setMethod("offTargets", "GuideSet", 
    function(object,
             columnName="alignments",
             max_mismatches=Inf,
             unlist=TRUE){
    if (!columnName %in% colnames(mcols(object)) ||
        !.isAlignmentsColumn(object, columnName)){
        out <- NULL
    } else {
        out <- mcols(object)[[columnName]]
        out <- BiocGenerics::unlist(out)
        out <- out[out$n_mismatches > 0]
        out <- out[out$n_mismatches <= max_mismatches]
        names(out) <- NULL
        if (!unlist){
            spacers <- spacers(object, as.character=TRUE)
            out <- split(out, f=factor(out$spacer,
                                       levels=unique(spacers)))
            out <- out[spacers]
            names(out) <- names(object)
        }
    }
    return(out)
})






#' @rdname GuideSet-class
#' @param gene_id Character vector of Ensembl gene IDs to subset gene
#'     annotation data by. If NULL (default), all genes are considered.
#' @param tx_id Character vector of Ensembl transcript IDs to subset gene
#'     annotation data by. If NULL (deafult), all transcript are considered.
#' @param gene_symbol Character vector of gene symbols to subset gene
#'     annotation data by. If NULL (default), all genes are considered.
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics unlist
#' @export
setMethod("geneAnnotation", "GuideSet", 
    function(object,
             unlist=TRUE,
             gene_id=NULL,
             tx_id=NULL,
             gene_symbol=NULL
             ){
    if (!"geneAnnotation" %in% colnames(S4Vectors::mcols(object))){
        out <- NULL
    } else {
        out <- S4Vectors::mcols(object)[["geneAnnotation"]]
        out <- BiocGenerics::unlist(out, use.names=FALSE)
        if (is.null(gene_id)){
            gene_id <- unique(out$gene_id)
        }
        if (is.null(tx_id)){
            tx_id <- unique(out$tx_id)
        }
        if (is.null(gene_symbol)){
            gene_symbol <- unique(out$gene_symbol)
        }
        cols <- c("gene_id", "tx_id", "gene_symbol")
        whs <- lapply(cols, function(col){
            out[[col]] %in% get(col)
        })
        wh <- Reduce("&", whs)
        out <- out[wh,,drop=FALSE]

        if (!unlist){
            out <- split(out, f=factor(rownames(out),
                                       levels=unique(names(object))))
            out <- out[names(object)]
            names(out) <- names(object)
        }
    }
    return(out)
})





#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics unlist
#' @export
setMethod("tssAnnotation", "GuideSet", 
    function(object,
             unlist=TRUE,
             gene_id=NULL,
             gene_symbol=NULL
             ){
    if (!"tssAnnotation" %in% colnames(S4Vectors::mcols(object))){
        out <- NULL
    } else {
        out <- S4Vectors::mcols(object)[["tssAnnotation"]]
        out <- BiocGenerics::unlist(out, use.names=FALSE)
        if (is.null(gene_id)){
            gene_id <- unique(out$gene_id)
        }

        if (is.null(gene_symbol)){
            gene_symbol <- unique(out$gene_symbol)     
        }
        cols <- c("gene_id", "gene_symbol")
        whs <- lapply(cols, function(col){
            out[[col]] %in% get(col)
        })
        wh <- Reduce("&", whs)
        out <- out[wh,,drop=FALSE]

        if (!unlist){
            out <- split(out, f=factor(rownames(out),
                                       levels=unique(names(object))))
            out <- out[names(object)]
            names(out) <- names(object)
        }
    }
    return(out)
})



#' @rdname GuideSet-class
#' @export
setMethod("enzymeAnnotation", "GuideSet", 
    function(object,
             unlist=TRUE){
    if (!"enzymeAnnotation" %in% colnames(mcols(object))){
        out <- NULL
        message("An enzymeAnnotation is not added yet. See ",
                "the function addRestrictionEnzymes to add ",
                "an enzyme annotation")
    } else {
        out <- mcols(object)[["enzymeAnnotation"]]
        if (unlist){
            out <- unlist(out)
        }
    }
    return(out)
})
