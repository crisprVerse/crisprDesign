#' @description Return \linkS4class{CrisprNuclease} object.
#' @rdname GuideSet-class
#' @export
setGeneric("crisprNuclease",
           function(object, ...) standardGeneric("crisprNuclease"))


#' @description Return string speecifying origin of DNA target.
#'     Either 'bsgenome' or 'customSequence'.
#' @rdname GuideSet-class
#' @export
setGeneric("targetOrigin",
           function(object, ...) standardGeneric("targetOrigin"))


#' @description Return custom DNA sequences used for designing gRNAs.
#' @rdname GuideSet-class
#' @export
setGeneric("customSequences",
           function(object, ...) standardGeneric("customSequences"))


#' @description Return BSgenome object used for designing gRNAs.
#' @rdname GuideSet-class
#' @export
setGeneric("bsgenome",
           function(object, ...) standardGeneric("bsgenome"))




#' @description Return spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("spacers",
           function(object, ...) standardGeneric("spacers"))


#' @description Return protospacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("protospacers",
           function(object, ...) standardGeneric("protospacers"))



#' @description Return PAM site coordinates.
#' @rdname GuideSet-class
#' @export
setGeneric("pamSites",
           function(object, ...) standardGeneric("pamSites"))


#' @description Return SNP annotation.
#' @rdname GuideSet-class
#' @export
setGeneric("snps",
           function(object, ...) standardGeneric("snps"))


#' @description Return genomic alignments of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("alignments",
           function(object, ...) standardGeneric("alignments"))



#' @description Return on-target alignments of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("onTargets",
           function(object, ...) standardGeneric("onTargets"))


#' @description Return off-target alignments of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("offTargets",
           function(object, ...) standardGeneric("offTargets"))



#' @description Return gene annotation table of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("geneAnnotation",
           function(object, ...) standardGeneric("geneAnnotation"))



#' @description Return TSS annotation table of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("tssAnnotation",
           function(object, ...) standardGeneric("tssAnnotation"))


#' @description Return restriction enzymes annotation table.
#' @rdname GuideSet-class
#' @export
setGeneric("enzymeAnnotation",
           function(object, ...) standardGeneric("enzymeAnnotation"))



