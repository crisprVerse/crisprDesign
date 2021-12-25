#' @title Add SNP annotation to a \linkS4class{GuideSet} object
#' 
#' @description Add SNP annotation to a \linkS4class{GuideSet} object.
#'    Only available for gRNAs designed for human genome.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param maf Minimum minor allele frequency to report
#'     (for a least one source among 1000Genomes and TOPMED).
#'     Default value of 0.01.
#' @param vcf Either string specfying path to a VCF file,
#'     or a \linkS4class{VCF} object.
#' @param verbose Should messages be printed? 
#'     TRUE by default.
#' 
#' @return The original \linkS4class{GuideSet} object with an additional
#'     list-column stored in \code{mcols{guideSet}} named \code{snpAnnotation}.
#' 
#' @seealso \code{link{snpAnnotation}} to retrieve an existing SNP annotation
#'     stored in a \linkS4class{GuideSet} object. See details section for a 
#'     description of the different columns.
#' 
#' @details
#' 
#' The different columns stored in
#' \code{mcols(guideSet)[["snps"]]} are:
#' 
#' \itemize{
#' \item \code{ID} gRNA ID.
#' \item \code{rs} Reference SNP cluster ID (e.g. rs17852242) 
#' \item \code{rs_site} Genomic coordinate of the SNP.
#' \item \code{rs_site_rel} SNP position relative to the PAM site.
#' \item \code{allele_ref} String specifying the reference allele of the SNP.
#' \item \code{allele_minor} String specifying the minor allele of the SNP.
#' \item \code{MAF_1000G} Minor allele frequency in the 1000 Genomes project.
#' \item \code{MAF_TOPMED} Minor allele frequency in the TOPMed project.
#' \item \code{snp_type} Character specifying the type of the SNP.
#' \item \code{snp_length} Character specifying the length of the SNP.
#' }
#' 
#' 
#' @examples
#'
#' data(guideSetExample, package="crisprDesign")
#' vcf <- system.file("extdata",
#'                    file="common_snps_dbsnp151_example.vcf.gz",
#'                    package="crisprDesign")
#' guideSet <- addSNPAnnotation(guideSetExample, vcf=vcf)
#' 
#' @export
#' @importFrom S4Vectors split
addSNPAnnotation <- function(guideSet,
                             maf=0.01,
                             vcf=NULL,
                             verbose=TRUE
){
    guideSet <- .validateGuideSet(guideSet)
    snps <- .getSNPAnnotation(guideSet=guideSet,
                              maf=maf,
                              vcf=vcf,
                              verbose=verbose)
    snps <- S4Vectors::split(snps,
                             f=factor(snps$ID,
                                      levels=names(guideSet)))
    mcols(guideSet)[["hasSNP"]] <- vapply(snps, nrow, FUN.VALUE=0)>0
    mcols(guideSet)[["snps"]] <- snps
    return(guideSet)
}




#' @importFrom VariantAnnotation ScanVcfParam readVcf info 
#' @importFrom MatrixGenerics rowRanges
#' @importFrom dplyr group_by bind_rows
#' @importFrom BiocGenerics strand
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors metadata DataFrame
.getSNPAnnotation <- function(guideSet,
                              maf=0.01,
                              vcf=NULL,
                              verbose=TRUE
){
    guideSet <- .validateGuideSet(guideSet)
    crisprNuclease <- crisprNuclease(guideSet)
    spacer_len <- spacerLength(crisprNuclease)
    genome <- metadata(guideSet)$genome
    if (genome!="hg38"){
        stop("SNPs can only be obtained for a GRanges in hg38 (human).")
    }
    vcf <- .validateVcf(vcf, verbose=verbose)

    # Checking maf argument
    stopifnot(is.numeric(maf))
    stopifnot(length(maf)==1)
    if (maf < 0 || maf > 1){
        stop('maf must be between 0 and 1.')
    }

    # To avoid chrM:
    valid_chrs <- paste0("chr", c(seq_len(22),"X","Y"))
    guideSet <- guideSet[as.character(seqnames(guideSet)) %in% valid_chrs]
    if (length(guideSet)==0){
        return(DataFrame())
    }

    # Prepare protospacer GRanges:
    gr.proto <- convertToProtospacerGRanges(guideSet)
    mcols(gr.proto) <- NULL
    seqnames <- seqnames(gr.proto)
    gr.proto <- .toggleSeqlevels(gr.proto)

    if (!is(vcf, "VCF")){
        svp <- ScanVcfParam(info=c('RS', 'RSPOS', 'CAF', 'TOPMED'),
                            which=gr.proto)
        vcf <- readVcf(vcf, genome='hg38',param=svp) 
    }

    # Get RS, RSPOS, CAF (1000G) and TOPMED:
    info <- info(vcf)
    rr   <- rowRanges(vcf)
    if (nrow(info)==0 || length(rr)==0){
        return(DataFrame())
    }
    ID <- as.character(rr$paramRangeID) # Spacer ID 
    rs <- paste0('rs', info$RS)         # SNP ID and position
    rs_site <- info$RSPOS               # SNP position

    # SNP position relative to PAM site:
    wh <- match(ID, names(guideSet))
    rs_site_rel   <- rs_site - guideSet$pam_site[wh]
    spacer_strand <- as.character(strand(guideSet)[wh])
    rs_site_rel[spacer_strand=='-'] <- -1 * rs_site_rel[spacer_strand=='-']

    # Reference and alt alleles:
    allele_ref   <- as.character(rr$REF)
    allele_minor <- vapply(rr$ALT, function(x){
        as.character(x[1])
    }, FUN.VALUE="a") 

    # Minor allele frequencies from 1000Genomes and TOPMED
    .getMinorFreq <- function(freq){
        freq <- vapply(freq, function(x){
            temp <- x[x!='.']
            temp <- as.numeric(temp)
            temp <- round(temp,5)
            temp[2]
        }, FUN.VALUE=1) 
        return(freq)
    }
    MAF.1000G  <- .getMinorFreq(info$CAF)
    MAF.TOPMED <- .getMinorFreq(info$TOPMED)

    # SNP type and length
    snp_type <- rep('snp', nrow(info))
    snp_type[nchar(allele_ref) > nchar(allele_minor)] <- 'del'
    snp_type[nchar(allele_ref) < nchar(allele_minor)] <- 'ins'
    snp_length <- nchar(allele_minor)

    # Get deletion size and insertion size:
    wh <- snp_type=='del'
    snp_length[wh] <- nchar(allele_ref[wh]) - nchar(allele_minor[wh])
    wh <- snp_type=='ins'
    snp_length[wh] <- nchar(allele_minor[wh])-nchar(allele_ref[snp_type=='ins'])

    # Add results to input data and return
    genome  <- "hg38"
    results <- DataFrame(ID = ID,
                         rs           = rs,
                         rs_site      = rs_site,
                         rs_site_rel  = rs_site_rel,
                         allele_ref   = allele_ref,
                         allele_minor = allele_minor,
                         MAF_1000G    = MAF.1000G,
                         MAF_TOPMED   = MAF.TOPMED,
                         snp_type     = snp_type,
                         snp_length   = snp_length)
    metadata(results)$crisprNuclease <- crisprNuclease
    metadata(results)$genome <- genome
    # Filtering on MAF:
    cond1 <- (!is.na(MAF.1000G)  & MAF.1000G  >= maf)
    cond2 <- (!is.na(MAF.TOPMED) & MAF.TOPMED >= maf)
    results <- results[cond1 | cond2,, drop=FALSE]
    return(results)
}






.validateVcf <- function(vcf, verbose=TRUE){
    if (is.character(vcf)){
        if (!file.exists(vcf)){
            msg <- paste0("Could not find the file ", vcf)
            stop(msg)
        }
    } else if (!is(vcf, "VCF")){
        stop("vcf must be either a VCF object or a path to a vcf file.")
    }
    return(vcf)
}




