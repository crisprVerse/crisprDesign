#library(crisprDesign)
#library(GenomicRanges)
#library(GenomeInfoDb)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
genome=BSgenome.Hsapiens.UCSC.hg38
seq <- getSeq(genome, "chr12", start=1, end=171330, as.character=TRUE)
seq <- DNAStringSet(seq)
names(seq) <- "chr12"
writeXStringSet(seq, "../fasta/chr12.fa")


#' # Exon-intro region of human KRAS:
#gr_input <- GRanges(c("chr12"),
#                    IRanges(start=25224014, end=25227007))
#genome(gr_input) <- "hg38"
#guideSet <- getGuides(gr_input)
#index="/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
#spacers <- spacers(guideSet, as.character=TRUE)
#aln <- getSpacerAlignments(spacers,
#                           bowtie_index=index,
#                           n_mismatches=1)
#aln <- aln[seqnames(aln) %in% paste0("chr", c(1:22, "X", "Y"))]
#aln <- GenomeInfoDb::keepStandardChromosomes(aln)
#aln <- promoters(aln, upstream=22, downstream=22)

#seqs <- getSeq(metadata(guideSet)$bsgenome,
#               aln,
#               as.character=TRUE)
#seq <- paste0(seqs, collapse="")
#seq <- DNAStringSet(seq)
#names(seq) <- "custom_chr"
#writeXStringSet(seq, "../fasta/custom_chr.fa")



