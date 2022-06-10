#' An S4 class to store pairs of CRISPR gRNA sequences.
#' 
#' @section Constructors:
#'     Use the constructor \code{link{PairedGuideSet}} to create a
#'     PairedGuideSet object.
#' 
#' @export
setClass("PairedGuideSet", contains="Pairs")


#' @importFrom S4Vectors Pairs
#' @export
PairedGuideSet <- function(GuideSet1 = NULL,
                           GuideSet2 = NULL
){
    .checkCutSitesOrder(GuideSet1, GuideSet2)
    pgs <- Pairs(GuideSet1, GuideSet2)
    mcols(pgs)[["pamOrientation"]] <- .getPamOrientation(GuideSet1,
                                                         GuideSet2)
    mcols(pgs)[["pamDistance"]] <- .getPamDistance(GuideSet1,
                                                   GuideSet2)
    mcols(pgs)[["spacerDistance"]] <- .getSpacerDistance(GuideSet1,
                                                         GuideSet2)
    mcols(pgs)[["cutLength"]] <- .getCutLength(GuideSet1,
                                               GuideSet2)
    new("PairedGuideSet", pgs)
}




setValidity("PairedGuideSet", function(object){
    out <- TRUE
    df <- mcols(object)
    mandatoryCols <- c("pamOrientation",
                       "pamDistance",
                       "spacerDistance",
                       "cutLength")
    
    if (!all(mandatoryCols %in% colnames(df))){
        diff <- setdiff(mandatoryCols, colnames(df))
        choices <- paste(diff, collapse=", ")
        out <- paste0("The following columns must be present",
               " in mcols(object): ",choices,".")
        return(out)
    } 

    # Checking pam orientation
    choices <- c("out", "in", "rev", "fwd")
    if (!all(df[["pamOrientation"]] %in% choices)){
        choices <- paste(choices, collapse=", ")
        out <- paste0("Only the following values are accepted ",
                      " pamOrientation: ",choices,".")
        return(out)
    }
    

    return(out)
})



### Check it's sorted correctly
.checkCutSitesOrder <- function(gs1,gs2){
    sites1 <- cutSites(gs1)
    sites2 <- cutSites(gs2)
    if (!all(sites1<sites2)){
        stop("Some of the pairs are not ordered properly. The cut site","
              of gRNA1 must be upstream of the cut site of gRNA2.")
    }
    return(TRUE)
}


#' @importFrom BiocGenerics strand
.getPamOrientation <- function(gs1, gs2){
    out <- rep(NA, length(gs1))
    s1 <- as.character(strand(gs1))
    s2 <- as.character(strand(gs2))
    out[s1=="-" & s2=="-"] <- "rev"
    out[s1=="+" & s2=="+"] <- "fwd"
    out[s1=="-" & s2=="+"] <- "out"
    out[s1=="+" & s2=="-"] <- "in"
    names(out) <- NULL
    return(out)
}



#' @importFrom crisprBase getProtospacerRanges
.getSpacerDistance <- function(gs1, gs2){
    r1 <- getProtospacerRanges(gs1,
                               nuclease=crisprNuclease(gs1))
    r2 <- getProtospacerRanges(gs2,
                               nuclease=crisprNuclease(gs2))
    out <- BiocGenerics::start(r2)-BiocGenerics::end(r1)
    out[!.onSameChr(gs1,gs2)] <- NA
    names(out) <- NULL
    return(out)
}


#' @importFrom crisprBase getProtospacerRanges
.getPamDistance <- function(gs1, gs2){
    r1 <- pamSites(gs1)
    r2 <- pamSites(gs2)
    out <- r2-r1
    names(out) <- NULL
    out[!.onSameChr(gs1,gs2)] <- NA
    return(out)
}




#' @importFrom crisprBase getProtospacerRanges
.getCutLength <- function(gs1, gs2){
    sites1 <- cutSites(gs1)
    sites2 <- cutSites(gs2)
    out <- sites2-sites1
    names(out) <- NULL
    out[!.onSameChr(gs1,gs2)] <- NA
    return(out)
}


#' @importFrom GenomeInfoDb seqnames
.onSameChr <- function(gs1, gs2){
    chr1 <- as.character(seqnames(gs1))
    chr2 <- as.character(seqnames(gs2))
    chr1==chr2
}



#' @rdname PairedGuideSet-class
#' @param ... Additional arguments
setMethod("pamOrientation", "PairedGuideSet",
    function(object){
    out <- mcols(object)[["pamOrientation"]]
    return(out)
})


#' @rdname PairedGuideSet-class
setMethod("pamDistance", "PairedGuideSet",
    function(object){
    out <- mcols(object)[["pamDistance"]]
    return(out)
})



#' @rdname PairedGuideSet-class
setMethod("pamDistance", "PairedGuideSet",
    function(object){
    out <- mcols(object)[["pamDistance"]]
    return(out)
})


#' @rdname PairedGuideSet-class
setMethod("cutLength", "PairedGuideSet",
    function(object){
    out <- mcols(object)[["cutLength"]]
    return(out)
})






#' @rdname PairedGuideSet-class
#' @param object \linkS4class{PairedGuideSet} object.
#' @param index Integer value indicating gRNA position.
#'     Must be either 1, 2, or NULL (default).
#'     If NULL, both positions are returned.
#' @export
#' @importFrom S4Vectors first second
#' @importFrom S4Vectors List
setMethod("crisprNuclease", "PairedGuideSet", 
    function(object, index=NULL){
    out <- List()
    out[[1]] <- metadata(first(object))[["CrisprNuclease"]]
    out[[2]] <- metadata(second(object))[["CrisprNuclease"]]
    names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    NULL
})




#' @rdname PairedGuideSet-class
#' @param as.character Should sequences be returned as a character
#'     vector? FALSE by default, in which case sequences are returned
#'     as a \linkS4class{DNAStringSet}.
#' @param returnAsRna Should the sequences be returned as RNA
#'     instead of DNA? FALSE by default. 
#' @export
#' @importFrom S4Vectors DataFrame
setMethod("spacers", "PairedGuideSet", 
    function(object,
             as.character=FALSE,
             returnAsRna=FALSE, 
             index=NULL){
    spacers1 <- spacers(first(object),
                        as.character=as.character,
                        returnAsRna=returnAsRna)
    spacers2 <- spacers(second(object),
                        as.character=as.character,
                        returnAsRna=returnAsRna)
    out <- DataFrame(first=spacers1,
                     second=spacers2)
    rownames(out) <- rownames(object)
    #names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    NULL
})




#' @rdname PairedGuideSet-class
#' @export
setMethod("pams", "PairedGuideSet", 
    function(object,
             as.character=FALSE,
             returnAsRna=FALSE,
             index=NULL){
    pams1 <- pams(first(object),
                  as.character=as.character,
                  returnAsRna=returnAsRna)
    pams2 <- pams(second(object),
                  as.character=as.character,
                  returnAsRna=returnAsRna)
    out <- DataFrame(first=pams1,
                     second=pams2)
    rownames(out) <- rownames(object)
    #names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    NULL
})




#' @rdname PairedGuideSet-class
#' @export
setMethod("pamSites", "PairedGuideSet", 
    function(object, index=NULL){
    sites1 <- pamSites(first(object))
    sites2 <- pamSites(second(object))
    out <- DataFrame(first=sites1,
                     second=sites2)
    rownames(out) <- rownames(object)
    #names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    return(out)
})




#' @rdname PairedGuideSet-class
#' @export
setMethod("cutSites", "PairedGuideSet", 
    function(object, index=NULL){
    sites1 <- cutSites(first(object))
    sites2 <- cutSites(second(object))
    out <- DataFrame(first=sites1,
                     second=sites2)
    rownames(out) <- rownames(object)
    #names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    return(out)
})




#' @rdname PairedGuideSet-class
#' @param include.pam Should PAM sequences be included?
#'     FALSE by default. 
#' @export
setMethod("protospacers", "PairedGuideSet", 
    function(object,
             as.character=FALSE,
             include.pam=FALSE,
             returnAsRna=FALSE,
             index=NULL){
    protospacers1 <- protospacers(first(object),
                                  as.character=as.character,
                                  include.pam=include.pam,
                                  returnAsRna=returnAsRna)
    protospacers2 <- protospacers(second(object),
                                  as.character=as.character,
                                  include.pam=include.pam,
                                  returnAsRna=returnAsRna)
    out <- DataFrame(first=protospacers1,
                     second=protospacers2)
    rownames(out) <- rownames(object)
    #names(out) <- c("first", "second")
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    NULL
})




#' @rdname PairedGuideSet-class
#' @export
setMethod("spacerLength", "PairedGuideSet", 
    function(object, index=NULL){
    nucs <- crisprNuclease(object)
    out <- unlist(lapply(nucs, spacerLength))
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    return(out)
})





#' @rdname PairedGuideSet-class
#' @export
setMethod("pamLength", "PairedGuideSet", 
    function(object, index=NULL){
    nucs <- crisprNuclease(object)
    out <- unlist(lapply(nucs, pamLength))
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    return(out)
})




#' @rdname PairedGuideSet-class
#' @export
setMethod("pamSide", "PairedGuideSet", 
    function(object, index=NULL){
    nucs <- crisprNuclease(object)
    out <- unlist(lapply(nucs, pamSide))
    if (is.null(index)){
        return(out)    
    } else if (index==1){
        return(out[[1]])
    } else if (index==2){
        return(out[[2]])
    } else {
        stop("index must be either NULL, 1, or 2.")
    }
    return(out)
})



# This merges the two GuideSet objects from a PairedGuideSet object
# into a unified GuideSet. This is useful to perform time intensive
# operations when there are a lot of duplicated guides between
# the two GuideSets.
.pairedGuideSet2GuideSet <- function(pairedGuideSet){
    nucs <- crisprNuclease(pairedGuideSet)
    nucsOK <- .identicalNucleases(nucs[[1]],
                                  nucs[[2]],
                                  checkWeights=TRUE,
                                  checkSpacerLength=TRUE)
    if (!nucsOK){
        stop("Cannot perform this operation as the gRNAs within pairs ",
             "do not originate from the same crisprNuclease object.")
    }
    gs1 <- .addCoordID(first(pairedGuideSet))
    gs2 <- .addCoordID(second(pairedGuideSet))
    gs <- unique(c(gs1,gs2))
    rownames(gs) <- NULL
    return(gs)
}


# This allows to add annotations to a PairedGuideSet object
# using annotations added to a unified GuideSet object
# obtained from the function .pairedGuideSet2GuideSet
#' @importFrom S4Vectors first<- second<-
.addColumnsFromUnifiedGuideSet <- function(pairedGuideSet,
                                           unifiedGuideSet
){
    first(pairedGuideSet) <- .borrowAnnotations(first(pairedGuideSet),
                                                unifiedGuideSet)
    second(pairedGuideSet) <- .borrowAnnotations(second(pairedGuideSet),
                                                 unifiedGuideSet)
    return(pairedGuideSet)
}

# Doing the heavy lifting for .addColumnsFromUnifiedGuideSet
.borrowAnnotations <- function(targetGuideSet,
                               sourceGuideSet
){
    nuc1 <- crisprNuclease(sourceGuideSet)
    nuc2 <- crisprNuclease(targetGuideSet)
    nucsOK <- .identicalNucleases(nuc1, nuc2,
                                  checkWeights=TRUE,
                                  checkSpacerLength=TRUE)
    if (!nucsOK){
        stop("Cannot perform this operation as the gRNAs within pairs ",
             "do not originate from the same crisprNuclease object.")
    }
    targetGuideSet <- .addCoordID(targetGuideSet)
    sourceGuideSet <- .addCoordID(sourceGuideSet)
    wh <- match(mcols(targetGuideSet)[["coordID"]],
                mcols(sourceGuideSet)[["coordID"]])
    if (sum(is.na(wh))>0){
        stop("Some gRNAs are not found in the sourceGuideSet.")
    }
    sourceGuideSet <- sourceGuideSet[wh]
    targetGuideSet_cols <- colnames(mcols(targetGuideSet))
    sourceGuideSet_cols <- colnames(mcols(sourceGuideSet))
    newCols <- setdiff(sourceGuideSet_cols,
                       targetGuideSet_cols)
    if (length(newCols)>0){
        for (k in seq_along(newCols)){
            col <- newCols[[k]]
            mcols(targetGuideSet)[[col]] <- mcols(sourceGuideSet)[[col]]
        }
    }
    return(targetGuideSet)
}




# library(crisprDesign)
# gs <- guideSetExample
# gs <- gs[order(BiocGenerics::start(gs))]
# gs1 <- c(gs[1:10],gs[1:10])
# gs2 <- c(gs[5:14+3],gs[5:14+80])
# pgs <- PairedGuideSet(gs1, gs2)
# gs <- .pairedGuideSetToGuideSet(pgs)


# pairedGuideSet <- pgs
# unifiedGuideSet <- .pairedGuideSetToGuideSet(pgs)
# unifiedGuideSet <- addPamScores(unifiedGuideSet)
# pairedGuideSet <- .addColumnsFromUnifiedGuideSet(pgs,
#                                                  unifiedGuideSet)











