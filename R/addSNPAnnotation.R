#' @title Add SNP annotation to a \linkS4class{GuideSet} object
#' 
#' @description Add SNP annotation to a \linkS4class{GuideSet} object.
#'    Only available for sgRNAs designed for human genome.
#' 
#' @param guideSet A \linkS4class{GuideSet} object.
#' @param vcf Either a character string specfying path to a VCF file
#'     or a \linkS4class{VCF} object.
#' @param maf Minimum minor allele frequency to report (for a least one source
#'     among 1000Genomes and TOPMED). Must be between 0 and 1 (exclusive).
#' 
#' @return \code{guideSet} appended with \code{hasSNP} column and \code{snps}
#'     list-column, both stored in \code{mcols{guideSet}}.
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
#' \item \code{ID} sgRNA ID.
#' \item \code{rs} Reference SNP cluster ID (e.g. rs17852242) 
#' \item \code{rs_site} Genomic coordinate of the SNP.
#' \item \code{rs_site_rel} Position of SNP relative to the PAM site.
#' \item \code{allele_ref} DNAString specifying the SNP reference allele.
#' \item \code{allele_minor} DNAString specifying the SNP minor allele.
#' \item \code{MAF_1000G} Minor allele frequency in the 1000 Genomes project.
#' \item \code{MAF_TOPMED} Minor allele frequency in the TOPMed project.
#' \item \code{type} Type of SNP (\code{"ins"}: insertion, \code{"del"}:
#'     deletion).
#' \item \code{length} Length of SNP in nucleotides.
#' }
#' 
#' @examples
#'
#' vcf <- system.file("extdata",
#'                    file="common_snps_dbsnp151_example.vcf.gz",
#'                    package="crisprDesign")
#' data(guideSetExample, package="crisprDesign")
#' guideSet <- addSNPAnnotation(guideSetExample, vcf=vcf)
#' 
#' 
#' 
#' @export
#' @importFrom S4Vectors split mcols<-
addSNPAnnotation <- function(guideSet,
                             vcf,
                             maf=0.01
){
    guideSet <- .validateGuideSet(guideSet)
    vcf <- .validateVcf(vcf, guideSet)
    stopifnot("maf must be a numeric value between 0 and 1" = {
        is.vector(maf, mode="numeric") &&
            length(maf) == 1 &&
            maf >= 0 &&
            maf < 1
    })
    snps <- .getSNPAnnotation(guideSet=guideSet,
                              maf=maf,
                              vcf=vcf)
    splitFactor <- factor(snps$ID, levels=unique(names(guideSet)))
    snps <- S4Vectors::split(snps, f=splitFactor)
    S4Vectors::mcols(guideSet)[["hasSNP"]] <- vapply(snps, function(x){
        nrow(x) > 0
    }, FUN.VALUE=logical(1))
    S4Vectors::mcols(guideSet)[["snps"]] <- snps
    return(guideSet)
}




#' @importFrom S4Vectors metadata metadata<- DataFrame mcols
#' @importFrom VariantAnnotation info 
#' @importFrom MatrixGenerics rowRanges
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings DNAStringSet
.getSNPAnnotation <- function(guideSet,
                              maf,
                              vcf
){
    genome <- S4Vectors::metadata(guideSet)[["genome"]]
    stopifnot("SNPs are only applicable for hg38 (human) genome" = {
        genome == "hg38"
    })
    
    # MtDNA not in tabix index
    validChrs <- paste0("chr", c(seq_len(22),"X","Y"))
    validChrs <- as.character(GenomeInfoDb::seqnames(guideSet)) %in% validChrs
    guideSet <- guideSet[validChrs]

    info <- VariantAnnotation::info(vcf)
    rr <- MatrixGenerics::rowRanges(vcf)
    rr <- S4Vectors::mcols(rr)
    
    if (!"paramRangeID" %in% colnames(rr)){
        rr <- cbind("paramRangeID"=character(0), rr)
    }
    snps <- cbind(rr, info)
    colnames(snps) <- c("ID", "allele_ref", "allele_minor", "rs", "rs_site",
                        "MAF_1000G", "MAF_TOPMED")
    
    snps$ID <- as.character(snps$ID)
    
    if (length(snps$allele_minor) > 0){
        snps$allele_minor <- lapply(snps$allele_minor, `[`, 1)
        snps$allele_minor <- do.call(c, snps$allele_minor)
    } else if (is(snps$allele_minor, "List")){
        snps$allele_minor <- unlist(snps$allele_minor)
    }
    
    snps$rs <- paste0("rs", snps$rs, recycle0=TRUE)
    
    guideIndices <- match(snps$ID, names(guideSet))
    rs_site_rel <- snps$rs_site - guideSet$pam_site[guideIndices]
    negativeStrand <- BiocGenerics::strand(guideSet)[guideIndices] == "-"
    negativeStrand <- as.logical(negativeStrand)
    rs_site_rel[negativeStrand] <- -1*rs_site_rel[negativeStrand]
    snps$rs_site_rel <- rs_site_rel
    
    # Minor allele frequencies from 1000Genomes and TOPMED
    .getMinorAlleleFreq <- function(freqs){
        vapply(freqs, function(x){
            minorFreq <- x[x != '.']
            minorFreq <- as.numeric(minorFreq)
            minorFreq <- round(minorFreq, 5)
            minorFreq[2]
        }, FUN.VALUE=numeric(1)) 
    }
    snps$MAF_1000G  <- .getMinorAlleleFreq(snps$MAF_1000G)
    snps$MAF_TOPMED <- .getMinorAlleleFreq(snps$MAF_TOPMED)
    
    # SNP type and length
    snpType <- rep("snp", nrow(snps))
    referenceLength <- nchar(snps$allele_ref)
    minorLength <- nchar(snps$allele_minor)
    snpDel <- referenceLength > minorLength
    snpIns <- referenceLength < minorLength
    snpType[snpDel] <- "del"
    snpType[snpIns] <- "ins"
    snps$type <- snpType
    
    snpLength <- nchar(snps$allele_ref)
    snpLength[snpDel] <- referenceLength[snpDel] - minorLength[snpDel]
    snpLength[snpIns] <- minorLength[snpIns] - referenceLength[snpIns]
    snps$length <- nchar(snps$allele_minor)
    
    snps <- snps[, c("ID", "rs", "rs_site", "rs_site_rel", "allele_ref",
                     "allele_minor", "MAF_1000G", "MAF_TOPMED", "type",
                     "length")]
    
    metadata(snps)$crisprNuclease <- crisprNuclease(guideSet)
    metadata(snps)$genome <- genome
    
    validMaf <- !is.na(snps$MAF_1000G) &
        !is.na(snps$MAF_TOPMED) &
        (snps$MAF_1000G >= maf | snps$MAF_TOPMED >= maf)
    
    snps <- snps[validMaf, , drop=FALSE]
    snps <- unique(snps)
    rownames(snps) <- snps$ID
    return(snps)
}





#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle<- seqinfo genome
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom VariantAnnotation ScanVcfParam readVcf
.validateVcf <- function(vcf,
                         guideSet
){
    if (is.character(vcf)){
        stopifnot("vcf must be a single file path" = {
            length(vcf) == 1
        })
        if (!file.exists(vcf)){
            missingFileMessage <- paste0("Could not find the file '", vcf, "'")
            stop(missingFileMessage)
        } else {
            protoGR <- convertToProtospacerGRanges(guideSet)
            if (length(protoGR) == 0){
                chrs <- VariantAnnotation::scanVcf(vcf)
                chrs <- GenomeInfoDb::seqlevels(chrs[[1]]$rowRanges)
                protoGR <- GenomicRanges::GRanges(
                    seqnames=chrs,
                    ranges=IRanges::IRanges(start=0, width=1))
            }
            GenomeInfoDb::seqlevelsStyle(protoGR) <- "Ensembl"
            genome <- GenomeInfoDb::seqinfo(protoGR)
            genome <- GenomeInfoDb::genome(genome)
            infoFields <- VariantAnnotation::scanVcfHeader(vcf)
            infoFields <- VariantAnnotation::info(infoFields)
            infoFields <- intersect(c("RS", "RSPOS", "CAF", "TOPMED"),
                                    rownames(infoFields))
            svp <- VariantAnnotation::ScanVcfParam(fixed="ALT",
                                                   info=infoFields,
                                                   which=protoGR)
            vcf <- VariantAnnotation::readVcf(vcf,
                                              genome=genome,
                                              param=svp)
        }
    } else {
        stopifnot("vcf must be either a VCF object or path to a vcf file" = {
            is(vcf, "VCF")
        })
    }
    return(vcf)
}



# fix rownames
