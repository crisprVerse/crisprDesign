############# Getter methods #############

#' @rdname GuideSet-class
#' @export
setGeneric("crisprNuclease",
           function(object, ...) standardGeneric("crisprNuclease"))


#' @rdname GuideSet-class
#' @export
setGeneric("targetOrigin",
           function(object, ...) standardGeneric("targetOrigin"))


#' @rdname GuideSet-class
#' @export
setGeneric("customSequences",
           function(object, ...) standardGeneric("customSequences"))



#' @rdname GuideSet-class
#' @export
setGeneric("bsgenome",
           function(object, ...) standardGeneric("bsgenome"))


#' @rdname GuideSet-class
#' @export
setGeneric("spacers",
           function(object, ...) standardGeneric("spacers"))


#' @rdname GuideSet-class
#' @export
setGeneric("protospacers",
           function(object, ...) standardGeneric("protospacers"))



#' @rdname GuideSet-class
#' @export
setGeneric("pamSites",
           function(object, ...) standardGeneric("pamSites"))


#' @rdname GuideSet-class
#' @export
setGeneric("snps",
           function(object, ...) standardGeneric("snps"))


#' @rdname GuideSet-class
#' @export
setGeneric("alignments",
           function(object, ...) standardGeneric("alignments"))


#' @rdname GuideSet-class
#' @export
setGeneric("onTargets",
           function(object, ...) standardGeneric("onTargets"))


#' @rdname GuideSet-class
#' @export
setGeneric("offTargets",
           function(object, ...) standardGeneric("offTargets"))


#' @rdname GuideSet-class
#' @export
setGeneric("geneAnnotation",
           function(object, ...) standardGeneric("geneAnnotation"))


#' @rdname GuideSet-class
#' @export
setGeneric("tssAnnotation",
           function(object, ...) standardGeneric("tssAnnotation"))


#' @rdname GuideSet-class
#' @export
setGeneric("enzymeAnnotation",
           function(object, ...) standardGeneric("enzymeAnnotation"))


#' @rdname GuideSet-class
#' @export
setGeneric("editedAlleles",
           function(object, ...) standardGeneric("editedAlleles"))


#' @rdname GuideSet-class
#' @export
setGeneric("txTable",
           function(object, ...) standardGeneric("txTable"))


#' @rdname GuideSet-class
#' @export
setGeneric("exonTable",
           function(object, ...) standardGeneric("exonTable"))


#' @rdname PairedGuideSet-class
#' @export
setGeneric("pamOrientation",
           function(object, ...) standardGeneric("pamOrientation"))


#' @rdname PairedGuideSet-class
#' @export
setGeneric("pamDistance",
           function(object, ...) standardGeneric("pamDistance"))



#' @rdname PairedGuideSet-class
#' @export
setGeneric("spacerDistance",
           function(object, ...) standardGeneric("spacerDistance"))


#' @rdname PairedGuideSet-class
#' @export
setGeneric("cutLength",
           function(object, ...) standardGeneric("cutLength"))




############# Setter methods #############

#' @rdname GuideSet-class
#' @export
setGeneric("tssAnnotation<-",
           function(object, value) standardGeneric("tssAnnotation<-"))


#' @rdname GuideSet-class
#' @export
setGeneric("geneAnnotation<-",
           function(object, value) standardGeneric("geneAnnotation<-"))


#' @rdname GuideSet-class
#' @export
setGeneric("enzymeAnnotation<-",
           function(object, value) standardGeneric("enzymeAnnotation<-"))


#' @rdname GuideSet-class
#' @export
setGeneric("snps<-",
           function(object, value) standardGeneric("snps<-"))

#' @rdname GuideSet-class
#' @export
setGeneric("alignments<-",
           function(object, value) standardGeneric("alignments<-"))




############# Add methods #############

#' @rdname addOnTargetScores
#' @export
setGeneric("addOnTargetScores",
           function(object, ...) standardGeneric("addOnTargetScores"))


#' @rdname addOffTargetScores
#' @export
setGeneric("addOffTargetScores",
           function(object, ...) standardGeneric("addOffTargetScores"))


#' @rdname addPamScores
#' @export
setGeneric("addPamScores",
           function(object, ...) standardGeneric("addPamScores"))


#' @rdname addCrispraiScores
#' @export
setGeneric("addCrispraiScores",
           function(object, ...) standardGeneric("addCrispraiScores"))


#' @rdname addCompositeScores
#' @export
setGeneric("addCompositeScores",
           function(object, ...) standardGeneric("addCompositeScores"))


#' @rdname addConservationScores
#' @export
setGeneric("addConservationScores",
           function(object, ...) standardGeneric("addConservationScores"))


#' @rdname addSNPAnnotation
#' @export
setGeneric("addSNPAnnotation",
           function(object, ...) standardGeneric("addSNPAnnotation"))


#' @rdname addIsoformAnnotation
#' @export
setGeneric("addIsoformAnnotation",
           function(object, ...) standardGeneric("addIsoformAnnotation"))



#' @rdname addDistanceToTss
#' @export
setGeneric("addDistanceToTss",
           function(object, ...) standardGeneric("addDistanceToTss"))




#' @rdname addRestrictionEnzymes
#' @export
setGeneric("addRestrictionEnzymes",
           function(object, ...) standardGeneric("addRestrictionEnzymes"))


#' @rdname addSequenceFeatures
#' @export
setGeneric("addSequenceFeatures",
           function(object, ...) standardGeneric("addSequenceFeatures"))


#' @rdname addRepeats
#' @export
setGeneric("addRepeats",
           function(object, ...) standardGeneric("addRepeats"))


#' @rdname removeRepeats
#' @export
setGeneric("removeRepeats",
           function(object, ...) standardGeneric("removeRepeats"))


#' @rdname addGeneAnnotation
#' @export
setGeneric("addGeneAnnotation",
           function(object, ...) standardGeneric("addGeneAnnotation"))


#' @rdname addTssAnnotation
#' @export
setGeneric("addTssAnnotation",
           function(object, ...) standardGeneric("addTssAnnotation"))


#' @rdname addPfamDomains
#' @export
setGeneric("addPfamDomains",
           function(object, ...) standardGeneric("addPfamDomains"))


#' @rdname addEditingSites
#' @export
setGeneric("addEditingSites",
           function(object, ...) standardGeneric("addEditingSites"))




#' @rdname GuideSet-class
#' @export
setGeneric("addCutSites",
           function(object, ...) standardGeneric("addCutSites"))





#' @rdname addSpacerAlignments
#' @export
setGeneric("addSpacerAlignments",
           function(object, ...) standardGeneric("addSpacerAlignments"))


#' @rdname addSpacerAlignments
#' @export
setGeneric("addSpacerAlignmentsIterative",
           function(object, ...) standardGeneric("addSpacerAlignmentsIterative"))




#' @rdname addNtcs
#' @export
setGeneric("addNtcs",
           function(object, ...) standardGeneric("addNtcs"))

