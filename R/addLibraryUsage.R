#' @title To add library annotation to a \linkS4class{GuideSet} object 
#' @description Function to add annotation to a \linkS4class{GuideSet} object
#'    to describe whether or not a given spacer sequence is found in an 
#'    existing CRISPR gRNA library. \pkg{crisprAnnotation} package must be
#'    installed. 
#' @param guideSet A \linkS4class{GuideSet} object.
#' 
#' @return Adds a column for each possible library present in the
#'    \pkg{crisprAnnotation} package that indicated whether or not
#'    a given spacer sequence is fount in the given library.
#' 
#' @author Luke Hoberecht
#' 
#' @examples 
#' \dontrun{
#'     guideSet <- addLibraryUsage(guideSetExample)
#' }
#' @export
#' @importFrom S4Vectors mcols<- mcols
addLibraryUsage <- function(guideSet){
    guideSet <- .validateGuideSet(guideSet)
    crisprNuclease <- crisprNuclease(guideSet)
    genome     <- metadata(guideSet)$genome
    spacer_len <- spacerLength(crisprNuclease)
    spacer_col <- paste0("spacer_", spacer_len, "mer")
    if (genome=="hg38"){
        species <- "human"
    } else if (genome=="mm10"){
        species <- "mouse"
    } else {
        return(guideSet)
    }
  
    if (!requireNamespace("crisprAnnotation")){
        warning("crisprAnnotation package is not installed. Skipping. ")
    } else {
        data(SpCas9, package="crisprBase", envir=environment())
        if (!identical(crisprNuclease, SpCas9)){
            message('Only libraries utilizing the SpCas9 nuclease ',
                    'are currently supported.', immediate.=TRUE)
            return(guideSet)
        }
        seqs <- spacers(guideSet,
                        as.character=TRUE)
    
        libs <- crisprAnnotation::listCrisprLibraries()
        if (species=="human"){
            pattern <- paste('avana',
                             'brunello',
                             'transylvania',
                             'sonata\\.wholegenome',
                             'human\\.gne[34]\\.wholegenome',
                             sep='|')
        } else {
            pattern <- paste('mouse\\.gne[34]\\.wholegenome',
                             sep='|')
        }
        libs <- libs[grepl(pattern=pattern, libs)]

        # get guides from libraries
        guides <- lapply(libs, function(x){
            temp <- crisprAnnotation::getCrisprLibraryAnnotation(x)
            if (spacer_col %in% colnames(temp)){
                temp <- temp[[spacer_col]]  
            } else {
                temp <- NULL
            }
            return(temp)
        })
        names(guides) <- libs
        guides <- compact(guides)
        
        # find matches in input data, looping over libraries
        if (length(guides)>0){
            for (i in seq_along(guides)){
                col_name <- paste0('lib_', names(guides)[[i]])
                wh <- seqs %in% guides[[i]]
                if (sum(wh)>0){
                    mcols(guideSet)[[col_name]] <- wh
                }
            }
        }
    }
    return(guideSet)
}