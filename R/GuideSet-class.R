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
#'     \item{\code{editedAlleles}:}{To get edited alleles annotation.}
#' }
#' 
#' @export
setClass("GuideSet", contains = "GRanges")






#' @describeIn GuideSet Create a \linkS4class{GuideSet} object
#' @param ids Character vector of unique gRNA ids. The ids can be anything,
#'     as long as they are unique. 
#' @param protospacers Character vector of protospacers sequences.
#' @param pams Character vector of PAM sequences.
#' @param seqnames Character vector of chromosome names.
#' @param pam_site Integer vector of PAM site coordinates.
#' @param strand Character vector of gRNA strand.
#'    Only accepted values are "+" and "-".
#' @param CrisprNuclease \linkS4class{CrisprNuclease} object.
#' @param targetOrigin String specifying the origin of the DNA target.
#'     Must be either 'bsgenome' or 'customSequences'.
#' @param bsgenome \linkS4class{BSgenome} object or string specifying
#'     BSgenome package name. Must be specified when
#'     \code{targetOrigin} is set to "bsgenome".
#' @param customSequences \linkS4class{DNAStringSet} object. Must be specified
#'     when \code{targetOrigin} is set to "customSequences".
#' @param ... Additional arguments for class-specific methods
#' @param seqinfo A \linkS4class{Seqinfo} object containing informatioon
#'     about the set of genomic sequences present in the target genome.
#' @param seqlengths \code{NULL}, or an integer vector named with
#'     \code{levels(seqnames)} and containing the lengths (or NA) for
#'     each level in \code{levels(seqnames)}.
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
#' ids <- paste0("grna_", seq_along(protospacers))
#' gr <- GuideSet(ids=ids,
#'                protospacers=protospacers,
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
GuideSet <- function(ids = NA_character_,
                     protospacers = NA_character_,
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
    
    # Checking ids
    if (sum(duplicated(ids))>0){
        stop("Duplicated values for 'ids' are not allowed.")
    }
    if (length(ids)!=length(protospacers)){
        stop("'ids' must have the same length as 'protospacers'.")
    }
    
    
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
        metadata(gr)[["bsgenome"]] <- .bsgenome_pkgname(bsgenome)
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
    names(gr) <- ids
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
        .isBSGenome(BSgenome::getBSgenome(meta[["bsgenome"]]))
    } else if (target=="customSequences"){
        .isDNAStringSet(meta[["customSequences"]])
    }
    
    nuc <- meta[["CrisprNuclease"]]
    if (!is(nuc, "CrisprNuclease")){
        out <- "metadata(object)$CrisprNuclease must be a CrisprNuclease object"
        return(out)
    }
    if (is.null(names(object))){
        out <- "GuideSet object cannot have NULL names."
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
#' @importFrom BSgenome getBSgenome
#' @export
setMethod("bsgenome", "GuideSet",
          function(object){
              out <- metadata(object)[["bsgenome"]]
              if (!is.character(out)){
                  out <- .bsgenome_pkgname(out)
              }
              out <- BSgenome::getBSgenome(out)
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
#' @importFrom crisprBase getCutSiteFromPamSite
setMethod("cutSites", "GuideSet", 
          function(object){
              pamSites <- mcols(object)[["pam_site"]]
              nuc <- metadata(object)[["CrisprNuclease"]]
              strand <- as.character(strand(object))
              ambiguousStrand <- strand == "*"
              out <- rep(NA, length(object))
              out[!ambiguousStrand] <- getCutSiteFromPamSite(
                  pam_site=pamSites[!ambiguousStrand],
                  strand=strand[!ambiguousStrand],
                  nuclease=nuc)
              names(out) <- names(object)
              return(out)
          })


#' @rdname GuideSet-class
#' @param object \linkS4class{GuideSet} object.
#' @export
#' @importFrom crisprBase getCutSiteFromPamSite
#' @importFrom S4Vectors mcols<-
setMethod("addCutSites", "GuideSet", 
          function(object){
              
              nuclease <- crisprNuclease(object)
              strand <- as.character(BiocGenerics::strand(object))
              ambiguousStrand <- strand == "*"
              cutSite <- rep(NA, length(object))
              cutSite[!ambiguousStrand] <- getCutSiteFromPamSite(
                  pam_site=pamSites(object)[!ambiguousStrand],
                  strand=strand[!ambiguousStrand],
                  nuclease=nuclease)
              mcols(object)[["cut_site"]] <- cutSite
              return(object)
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
#' @param use.names Whether to include spacer IDs as (row)names (\code{TRUE}),
#'     or as a separate column (\code{FALSE}).
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("snps", "GuideSet", 
          function(object,
                   unlist=TRUE,
                   use.names=TRUE){
              if (!"snps" %in% colnames(S4Vectors::mcols(object))){
                  out <- NULL
              } else {
                  out <- S4Vectors::mcols(object)[["snps"]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          })


#' @rdname GuideSet-class
#' @param columnName Name of the column storing the alignments annotation
#'     to be retrieved.
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist
#' @export
setMethod("alignments", "GuideSet", 
          function(object,
                   columnName="alignments",
                   unlist=TRUE,
                   use.names=TRUE){
              if (!columnName %in% colnames(S4Vectors::mcols(object)) ||
                  !.isAlignmentsColumn(object, columnName)){
                  out <- NULL
              } else {
                  out <- S4Vectors::mcols(object)[[columnName]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!unlist){
                      split_factor <- factor(names(out), levels=names(object))
                      if (!use.names){
                          out <- .namesAsColumn_gr(out)
                      }   
                      out <- S4Vectors::split(out, f=split_factor)[names(object)]
                  }
              }
              return(out)
          })





#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist
#' @export
setMethod("onTargets", "GuideSet", 
          function(object,
                   columnName="alignments",
                   unlist=TRUE,
                   use.names=TRUE){
              if (!columnName %in% colnames(S4Vectors::mcols(object)) ||
                  !.isAlignmentsColumn(object, columnName)){
                  out <- NULL
              } else {
                  out <- S4Vectors::mcols(object)[[columnName]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  out <- out[out$n_mismatches == 0]
                  if (!unlist){
                      split_factor <- factor(names(out), levels=names(object))
                      if (!use.names){
                          out <- .namesAsColumn_gr(out)
                      } 
                      out <- S4Vectors::split(out, f=split_factor)[names(object)]
                  }
              }
              return(out)
          })


#' @rdname GuideSet-class
#' @param max_mismatches What should be the maximum number of 
#'     mismatches considered for off-targets? 
#'     Inf by default.
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist
#' @export
setMethod("offTargets", "GuideSet", 
          function(object,
                   columnName="alignments",
                   max_mismatches=Inf,
                   unlist=TRUE,
                   use.names=TRUE){
              
              stopifnot("max_mismatches must be a non-negative integer" = {
                  is.vector(max_mismatches, mode="numeric") &&
                      length(max_mismatches) == 1 &&
                      max_mismatches == round(max_mismatches) &&
                      max_mismatches >= 0
              })
              if (!columnName %in% colnames(S4Vectors::mcols(object)) ||
                  !.isAlignmentsColumn(object, columnName)){
                  out <- NULL
              } else {
                  out <- S4Vectors::mcols(object)[[columnName]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  out <- out[out$n_mismatches > 0 & out$n_mismatches <= max_mismatches]
                  if (!unlist){
                      split_factor <- factor(names(out), levels=names(object))
                      if (!use.names){
                          out <- .namesAsColumn_gr(out)
                      } 
                      out <- S4Vectors::split(out, f=split_factor)[names(object)]
                  }
              }
              return(out)
          })





# setMethod("alignments", "GuideSet", 
#     function(object,
#              columnName="alignments",
#              unlist=TRUE,
#              use.names=TRUE){
#     if (!columnName %in% colnames(S4Vectors::mcols(object)) ||
#         !.isAlignmentsColumn(object, columnName)){
#         out <- NULL
#     } else {
#         out <- S4Vectors::mcols(object)[[columnName]]
#         out <- BiocGenerics::unlist(out, use.names=FALSE)
#         if (!use.names){
#             out <- .namesAsColumn_gr(out)
#             split_factor <- S4Vectors::mcols(out)[["spacer_id"]]
#         } else {
#             split_factor <- names(out)
#         }
#         if (!unlist){
#             out <- S4Vectors::split(out, f=split_factor)
#         }
#     }
#     return(out)
# })




#' @rdname GuideSet-class
#' @param value Object to replace with
#' @export
setMethod("alignments<-", "GuideSet", 
          function(object, value){
              mcols(object)[["alignments"]] <- value
              return(object)
          })

#' @rdname GuideSet-class
#' @export
setMethod("geneAnnotation<-", "GuideSet", 
          function(object, value){
              mcols(object)[["geneAnnotation"]] <- value
              return(object)
          })

#' @rdname GuideSet-class
#' @export
setMethod("tssAnnotation<-", "GuideSet", 
          function(object, value){
              mcols(object)[["tssAnnotation"]] <- value
              return(object)
          })

#' @rdname GuideSet-class
#' @export
setMethod("enzymeAnnotation<-", "GuideSet", 
          function(object, value){
              mcols(object)[["enzymeAnnotation"]] <- value
              return(object)
          })


#' @rdname GuideSet-class
#' @export
setMethod("snps<-", "GuideSet", 
          function(object, value){
              mcols(object)[["snps"]] <- value
              return(object)
          })


# #' @rdname GuideSet-class
# #' @export
# setMethod("txTable<-", "GuideSet", 
#     function(object, value){
#     mcols(object)[["txTable"]] <- value
#     return(object)
# })












#' @rdname GuideSet-class
#' @param gene_id Character vector of Ensembl gene IDs to subset gene
#'     annotation data by. If NULL (default), all genes are considered.
#' @param tx_id Character vector of Ensembl transcript IDs to subset gene
#'     annotation data by. If NULL (deafult), all transcript are considered.
#' @param gene_symbol Character vector of gene symbols to subset gene
#'     annotation data by. If NULL (default), all genes are considered.
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("geneAnnotation", "GuideSet", 
          function(object,
                   unlist=TRUE,
                   gene_id=NULL,
                   tx_id=NULL,
                   gene_symbol=NULL,
                   use.names=TRUE){
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
                  cols <- intersect(cols, colnames(out))
                  whs <- lapply(cols, function(col){
                      out[[col]] %in% get(col)
                  })
                  wh <- Reduce("&", whs)
                  out <- out[wh, , drop=FALSE]
                  
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          })




#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols
#' @export
setMethod("editedAlleles", "GuideSet", 
          function(object,
                   unlist=TRUE,
                   use.names=TRUE){
              if (!"editedAlleles" %in% colnames(S4Vectors::mcols(object))){
                  out <- NULL
                  message("Edited alleles annotation has not been added yet.",
                          "See the function 'addEditedAlleles' to add ",
                          "edited alleles annotation.")
              } else {
                  out <- S4Vectors::mcols(object)[["editedAlleles"]]
                  out <- do.call(rbind, out)
                  # out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          })




#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("tssAnnotation", "GuideSet", 
          function(object,
                   unlist=TRUE,
                   gene_id=NULL,
                   gene_symbol=NULL,
                   use.names=TRUE){
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
                  out <- out[wh, , drop=FALSE]
                  
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          })



#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("enzymeAnnotation", "GuideSet", 
          function(object,
                   unlist=TRUE,
                   use.names=TRUE){
              if (!"enzymeAnnotation" %in% colnames(S4Vectors::mcols(object))){
                  out <- NULL
                  message("An enzymeAnnotation is not added yet. See ",
                          "the function 'addRestrictionEnzymes' to add ",
                          "enzyme annotation")
              } else {
                  out <- S4Vectors::mcols(object)[["enzymeAnnotation"]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          })




#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("txTable", "GuideSet",
          function(object,
                   unlist=TRUE,
                   use.names=TRUE){
              if (!"txTable" %in% colnames(S4Vectors::mcols(object))){
                  out <- NULL
                  message("A txTable has not been added yet. See ",
                          "the function 'addTxTable' to add ",
                          "a txTable.")
              } else {
                  out <- S4Vectors::mcols(object)[["txTable"]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          }
)




#' @rdname GuideSet-class
#' @importFrom S4Vectors mcols split
#' @importFrom BiocGenerics unlist rownames
#' @export
setMethod("exonTable", "GuideSet",
          function(object,
                   unlist=TRUE,
                   use.names=TRUE){
              if (!"exonTable" %in% colnames(S4Vectors::mcols(object))){
                  out <- NULL
                  message("An exonTable has not been added yet. See ",
                          "the function 'addExonTable' to add ",
                          "a txTable.")
              } else {
                  out <- S4Vectors::mcols(object)[["exonTable"]]
                  out <- BiocGenerics::unlist(out, use.names=FALSE)
                  if (!use.names){
                      out <- .namesAsColumn_df(out)
                      split_factor <- out[["spacer_id"]]
                  } else {
                      split_factor <- BiocGenerics::rownames(out)
                  }
                  if (!unlist){
                      out <- S4Vectors::split(out, f=split_factor)
                  }
              }
              return(out)
          }
)







# Add data.frame rownames to a column 
#' @importFrom BiocGenerics rownames rownames<- cbind
.namesAsColumn_df <- function(out){
    spacer_id <- BiocGenerics::rownames(out)
    out <- BiocGenerics::cbind(spacer_id, out)
    BiocGenerics::rownames(out) <- NULL
    return(out)
}


# Add GRanges names to a metadata column 
#' @importFrom S4Vectors DataFrame mcols mcols<-
.namesAsColumn_gr <- function(out
){
    spacer_id <- names(out)
    new_mcols <- S4Vectors::DataFrame(spacer_id,
                                      S4Vectors::mcols(out))
    S4Vectors::mcols(out) <- new_mcols
    names(out) <- NULL
    return(out)
}
