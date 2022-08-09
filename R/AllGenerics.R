############# Getter methods #############

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


#' @description Return list of edited alleles
#' @rdname GuideSet-class
#' @export
setGeneric("editedAlleles",
           function(object, ...) standardGeneric("editedAlleles"))

#' @description Return PAM orientation configuration
#' @rdname PairedGuideSet-class
#' @export
setGeneric("pamOrientation",
           function(object, ...) standardGeneric("pamOrientation"))


#' @description Return distance between PAM sites from paired gRNAs
#' @rdname PairedGuideSet-class
#' @export
setGeneric("pamDistance",
           function(object, ...) standardGeneric("pamDistance"))



#' @description Return distance between spacer sequences from paired gRNAs
#' @rdname PairedGuideSet-class
#' @export
setGeneric("spacerDistance",
           function(object, ...) standardGeneric("spacerDistance"))


#' @description Return cut length resulting from paired gRNAs
#' @rdname PairedGuideSet-class
#' @export
setGeneric("cutLength",
           function(object, ...) standardGeneric("cutLength"))




############# Setter methods #############

#' @description Store TSS annotation table of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("tssAnnotation<-",
           function(object, value) standardGeneric("tssAnnotation<-"))


#' @description Store gene annotation table of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("geneAnnotation<-",
           function(object, value) standardGeneric("geneAnnotation<-"))


#' @description Store restriction enzymes annotation table.
#' @rdname GuideSet-class
#' @export
setGeneric("enzymeAnnotation<-",
           function(object, value) standardGeneric("enzymeAnnotation<-"))


#' @description Store SNP annotation.
#' @rdname GuideSet-class
#' @export
setGeneric("snps<-",
           function(object, value) standardGeneric("snps<-"))

#' @description Store genomic alignments of spacer sequences.
#' @rdname GuideSet-class
#' @export
setGeneric("alignments<-",
           function(object, value) standardGeneric("alignments<-"))




############# Add methods #############

#' @description Add on-target scores.
#' @rdname GuideSet-class
#' @export
setGeneric("addOnTargetScores",
           function(object, ...) standardGeneric("addOnTargetScores"))


#' @description Add off-target scores.
#' @rdname GuideSet-class
#' @export
setGeneric("addOffTargetScores",
           function(object, ...) standardGeneric("addOffTargetScores"))




#' @description Add PAM scores.
#' @rdname GuideSet-class
#' @export
setGeneric("addPamScores",
           function(object, ...) standardGeneric("addPamScores"))


#' @description Add CRISPRai scores.
#' @rdname GuideSet-class
#' @export
setGeneric("addCrispraiScores",
           function(object, ...) standardGeneric("addCrispraiScores"))



#' @description Add composite scores.
#' @rdname GuideSet-class
#' @export
setGeneric("addCompositeScores",
           function(object, ...) standardGeneric("addCompositeScores"))





#' @description Add SNP annotation.
#' @rdname GuideSet-class
#' @export
setGeneric("addSNPAnnotation",
           function(object, ...) standardGeneric("addSNPAnnotation"))



#' @description Add restriction enzymes annotation.
#' @rdname GuideSet-class
#' @export
setGeneric("addRestrictionEnzymes",
           function(object, ...) standardGeneric("addRestrictionEnzymes"))


#' @description Add spacer sequence features.
#' @rdname GuideSet-class
#' @export
setGeneric("addSequenceFeatures",
           function(object, ...) standardGeneric("addSequenceFeatures"))

#' @description Add spacer sequence features.
#' @rdname GuideSet-class
#' @export
setGeneric("addSequenceFeatures",
           function(object, ...) standardGeneric("addSequenceFeatures"))


#' @description Annotate with repeat elements
#' @rdname GuideSet-class
#' @export
setGeneric("addRepeats",
           function(object, ...) standardGeneric("addRepeats"))


#' @description Remove rows overlaping repeat elements
#' @rdname GuideSet-class
#' @export
setGeneric("removeRepeats",
           function(object, ...) standardGeneric("removeRepeats"))


#' @description Add gene context annotation
#' @rdname GuideSet-class
#' @export
setGeneric("addGeneAnnotation",
           function(object, ...) standardGeneric("addGeneAnnotation"))


#' @description Add TSS context annotation
#' @rdname GuideSet-class
#' @export
setGeneric("addTssAnnotation",
           function(object, ...) standardGeneric("addTssAnnotation"))


#' @description Add spacer on- and off-target alignments
#' @rdname GuideSet-class
#' @export
setGeneric("addSpacerAlignments",
           function(object, ...) standardGeneric("addSpacerAlignments"))


#' @description Add spacer on- and off-target alignments (iterative mode)
#' @rdname GuideSet-class
#' @export
setGeneric("addSpacerAlignmentsIterative",
           function(object, ...) standardGeneric("addSpacerAlignmentsIterative"))


