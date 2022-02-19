#' @title To add edited alleles resulting from base editing
#' 
#' @description To add edited alleles resulting from base editing
#' 
#' @param txInfoDataFrame A \code{DataFrame} obtained from \code{getTxInfoDataFrame}.
#' @param substitution String specifying the base editing substitution.
#'     E.g: "C2T" to specify that Cs are converted to Ts.
#' 
#' @return The original \code{txInfoDataFrame} object with the following
#'     additional columns:
#' 
#' \itemize{
#' \item \code{aa_edited_fwd} Amino acid resulting from editing forward strand.
#' \item \code{nuc_edited_fwd}  Nucleotide resulting from editing forward strand.
#' \item \code{class_fwd} Type of mutation (silent, missense, or nonsense)
#'     resulting from editing forward strand,
#' \item \code{stop_fwd} Does editing on the forward strand result in a stop codon?
#' \item \code{aa_edited_rev} Amino acid resulting from editing reverse strand.
#' \item \code{nuc_edited_rev} Nucleotide resulting from editing reverse strand.
#' \item \code{class_rev} Type of mutation (silent, missense, or nonsense)
#'     resulting from editing forward strand.
#' \item\code{stop_rev} Does editing on the reverse strand result in a stop codon?
#' }
#' 
#' @examples 
#' 
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#' data("grListExample")
#' tx_id <- "ENST00000538872"
#' df <- getTxInfoDataFrame(tx_id=tx_id,
#'     txObject=grListExample,
#'     bsgenome=bsgenome)
#' df <- addEditedAlleles(df, substitution="C2T")
#' 
#' @author Jean-Philippe Fortin
#' 
#' @importFrom Biostrings DNAString translate
#' @export
addEditedAlleles <- function(txInfoDataFrame,
                             substitution="C2T"
){
    if (length(substitution)>1){
        stop("substitution must be a string of length 1")
    }
    substitution <- .validateSubstitution(substitution)


    gene_strand <- txInfoDataFrame$strand[[1]]
    cds <- txInfoDataFrame[!is.na(txInfoDataFrame$aa_number),,drop = FALSE]
    dna <- paste0(cds$nuc, collapse="")
    dna <- DNAString(dna)
    protein <- translate(dna)
    protein <- strsplit(as.character(protein), "")[[1]]
    protein <- rep(protein, each=3)
    
    # On original strand:
    .addEditing <- function(cds,
                            strand=c("original", "reverse")){
        strand <- match.arg(strand)    
        editingStrand <- .getEditingStrand(strand=strand,
                                           gene_strand=gene_strand)
        originBase <- .getOriginBase(substitution, strand)
        targetBase <- .getTargetBase(substitution, strand)
        dna_edited <- DNAString(gsub(originBase, targetBase, dna))
        protein_edited <- translate(dna_edited)
        protein_edited <- strsplit(as.character(protein_edited), "")[[1]]
        protein_edited <- rep(protein_edited, each=3)
        dna_edited <- strsplit(as.character(dna_edited), "")[[1]]
        if (editingStrand=="fwd"){
            cds$aa_edited_fwd  <- protein_edited
            cds$nuc_edited_fwd <- dna_edited
        } else {
            cds$aa_edited_rev  <- protein_edited
            cds$nuc_edited_rev <- dna_edited
        }
        return(cds)
    }
    cds <- .addEditing(cds, "original") 
    cds <- .addEditing(cds, "reverse") 
    cds$class_fwd <- NA
    cds$class_rev <- NA

    # Annotating mutations:
    cds$class_fwd[cds$aa_edited_fwd==cds$aa]  <- "silent"
    cds$class_rev[cds$aa_edited_rev==cds$aa]  <- "silent"
    cds$class_fwd[cds$aa_edited_fwd!=cds$aa]  <- "missense"
    cds$class_rev[cds$aa_edited_rev!=cds$aa]  <- "missense"
    cds$class_fwd[cds$aa_edited_fwd=="*"]     <- "nonsense"
    cds$class_rev[cds$aa_edited_rev=="*"]     <- "nonsense"

    cols <- c("pos_mrna", 
              "aa_edited_fwd",
              "nuc_edited_fwd",
              "class_fwd",
              "aa_edited_rev", 
              "nuc_edited_rev",
              "class_rev")
    cds <- cds[,cols, drop=FALSE]
    wh <- match(txInfoDataFrame$pos_mrna, cds$pos_mrna)
    cds <- cds[wh,,drop=FALSE]
    cds$pos_mrna <- NULL
    txInfoDataFrame <- cbind(txInfoDataFrame, cds)
    rownames(txInfoDataFrame) <- NULL

    # Adding stop codons:
    txInfoDataFrame$stop_fwd <- txInfoDataFrame$aa_edited_fwd=="*"
    txInfoDataFrame$stop_rev <- txInfoDataFrame$aa_edited_rev=="*"
    txInfoDataFrame$stop_fwd[is.na(txInfoDataFrame$stop_fwd)] <- FALSE
    txInfoDataFrame$stop_rev[is.na(txInfoDataFrame$stop_rev)] <- FALSE

    return(txInfoDataFrame)
}





.getBaseEditingKey <- function(){
    
    .getComboNames <- function(){
        dnaLetters <- c("A", "C", "G", "T")
        combs <- expand.grid(dnaLetters, dnaLetters)
        combs <- combs[combs[,1]!=combs[,2],]
        combs <- paste0(combs[,1], "2", combs[,2])
        return(combs)
    }


    key <- data.frame(substitution=.getComboNames())
    key$origin <- .getOriginBase(key$substitution)
    key$target <- .getTargetBase(key$substitution)
    return(key)
}




#' @importFrom Biostrings complement
.getOriginBase <- function(x,
                           strand=c("original","reverse")
){
    strand <- match.arg(strand)
    x <- substr(x,1,1)
    if (strand=="reverse"){
        x <- complement(DNAStringSet(x))
        x <- as.character(x)
    }
    return(x)
}

#' @importFrom Biostrings complement
.getTargetBase <- function(x,
                           strand=c("original","reverse")
){
    strand <- match.arg(strand)
    x <- substr(x,3,3)
    if (strand=="reverse"){
        x <- complement(DNAStringSet(x))
        x <- as.character(x)
    }
    return(x)
}


.validateSubstitution <- function(substitution){
    substitution <- toupper(substitution)
    beKey <- .getBaseEditingKey()
    choices <- beKey$substitution
    if (!substitution %in% choices){
        choices <- paste0(choices, collapse=", ")
        choices <- paste0(choices, ".")
        stop("The specified substitution is not valid.",
             "It must be one of the following: ",
             choices)
    }
    return(substitution)
}



.getEditingStrand <- function(strand,
                              gene_strand
){
    if (strand=="original" & gene_strand=="+"){
        out <- "fwd"
    } else if (strand=="original" & gene_strand=="-"){
        out <- "rev"
    } else if (strand=="reverse" & gene_strand=="-"){
        out <- "fwd"
    } else {
        out <- "rev"
    }
    return(out)
}





# #' @title To annotate stop codons resulting from base editing
# #' 
# #' @description To annotate stop codons resulting from base editing.
# #' 
# #' @param txInfoDataFrame A \code{DataFrame} obtained from \code{getTxInfoDataFrame}.
# #' @param substitution String specifying the base editing substitution.
# #'     E.g: "C2T" to specify that Cs are converted to Ts.
# #' 
# #' @return The original \code{txInfoDataFrame} object with the following
# #'     additional columns: \code{stop_fwd} and \code{stop_rev}.
# #'     \code{stop_fwd} is a logical vector indicating whether or not
# #'     a stop codon is created by editing the forward strand of 
# #'     the target DNA, and \code{stop_rev} is a logical value indicating
# #'     whether or not a stop codon is created by editing the reverse
# #'     strand of the target DNA. 
# #' 
# #'     Gene strand was taken into account when generating those columns.
# #'     For instance, \code{stop_fwd} equal to TRUE indicates that base editing
# #'     occuring on the forward strand of the DNA will esult in a stop
# #'     codon in the specified gene, whether or not the gene is located
# #'     on the forward or reverse strand. 
# #' 
# #' @examples 
# #' 
# #' library(BSgenome.Hsapiens.UCSC.hg38)
# #' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
# #' data("grListExample")
# #' tx_id <- "ENST00000538872"
# #' df <- getTxInfoDataFrame(tx_id=tx_id,
# #'     txObject=grListExample,
# #'     bsgenome=bsgenome)
# #' df <- annotateStopCodons(df, substitution="C2T")
# #' 
# #' @author Jean-Philippe Fortin
# #' 
# #' @export
# annotateStopCodons <- function(txInfoDataFrame,
#                                substitution="C2T"
# ){
#     if (length(substitution)>1){
#         stop("substitution must be a string of length 1")
#     }
#     substitution <- .validateSubstitution(substitution)
    
#     gene_strand <- txInfoDataFrame$strand[[1]]

#     # Only considering CDS:
#     cds <- txInfoDataFrame[!is.na(txInfoDataFrame$aa_number),,drop=FALSE]
#     cds$stop_fwd <- FALSE
#     cds$stop_rev <- FALSE
#     codons <- split(cds, f=cds$aa_number)

#     .annotateCodonsByStrand <- function(codons,
#                                         strand){
#         codons <- lapply(codons, function(x){
#             y <- x$nuc
#             originBase <- .getOriginBase(substitution, strand)
#             targetBase <- .getTargetBase(substitution, strand)
#             y[y==originBase] <- targetBase
            
#             editingStrand <- .getEditingStrand(strand=strand,
#                                                gene_strand=gene_strand)
#             codon <- paste0(y, collapse="")
#             if (codon %in% STOP_CODONS){
#                 if (editingStrand=="fwd"){
#                     x$stop_fwd <- TRUE
#                 } else {
#                     x$stop_rev <- TRUE
#                 }
#             } 
#             return(x)
#         })
#         return(codons)
#     }    
#     codons <- .annotateCodonsByStrand(codons, "original")
#     codons <- .annotateCodonsByStrand(codons, "reverse")
#     codons <- do.call("rbind", codons)
#     txInfoDataFrame$stop_fwd <- txInfoDataFrame$pos %in% codons$pos[codons$stop_fwd]
#     txInfoDataFrame$stop_rev <- txInfoDataFrame$pos %in% codons$pos[codons$stop_rev]
#     return(txInfoDataFrame)
# }





















