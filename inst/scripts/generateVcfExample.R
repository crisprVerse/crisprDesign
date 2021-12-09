library(VariantAnnotation)
library(GenomicRanges)
library(GenomeInfoDb)
inputfile <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"

load("../../data/grListExample.rda")
min <- min(start(grListExample$cds))
max <- max(end(grListExample$cds))
chr <- as.character(seqnames(grListExample$cds)[1])
chr <- gsub("chr","", chr)
gr <- GRanges(chr, IRanges(start=min, end=max))
genome(gr) <- "hg38"

#svp <- ScanVcfParam(info=c("CAF"))
svp <- ScanVcfParam(which=gr, info=c("CAF"))
svp <- ScanVcfParam(which=gr)
vcf <- readVcf(inputfile,
               param=svp,
               genome=genome(gr))
writeVcf(vcf,
         file="../extdata/common_snps_dbsnp151_example.vcf")

#bgzip common_snps_dbsnp151_example.vcf
#tabix -p vcf common_snps_dbsnp151_example.vcf.gz