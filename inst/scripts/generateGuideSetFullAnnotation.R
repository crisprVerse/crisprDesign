library(crisprDesign)
library(crisprDesignData)
library(crisprScore)
library(biomaRt)
library(devtools)
library(BSgenome.Hsapiens.UCSC.hg38)
data(tss_human, package="crisprDesignData")

txObject  <- txdb_human
tssObject <- tss_human
gr.repeats <- gr.repeats.hg38
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
# bowtie_index <- "~/crisprIndices/bowtie/hg38/hg38"
bwa_index <- "~/crisprIndices/bwa/hg38/hg38"
vcf <- "~/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"
mart_dataset <- "hsapiens_gene_ensembl"

load("../../data/guideSetExample.rda")
gs <- guideSetExample
set.seed(46076919)
indices <- sample(length(gs), 20)
gs <- gs[indices]


gs <- addSequenceFeatures(gs,
                          addHairpin=TRUE)
gs <- addSpacerAlignmentsIterative(gs,
                                   aligner='bwa',
                                   txObject=txObject,
                                   tssObject=tssObject,
                                   n_mismatches=3,
                                   aligner_index=bwa_index,
                                   bsgenome=bsgenome)
guideSetExampleWithAlignments <- gs
use_data(guideSetExampleWithAlignments,
         compress="xz",
         overwrite=TRUE)


gs <- addOffTargetScores(gs)
gs <- addOnTargetScores(gs)
gs <- addPamScores(gs)
gs <- addRepeats(gs,
                 gr.repeats=gr.repeats.hg38)
gs <- addRestrictionEnzymes(gs)
gs <- addGeneAnnotation(gs,
                        txObject=txObject,
                        ignore_introns=FALSE,
                        addPfam=TRUE,
                        mart_dataset=mart_dataset)
gs <- addTssAnnotation(gs,
                       tssObject=tssObject)
gs <- addSNPAnnotation(gs, vcf=vcf)



guideSetExampleFullAnnotation <- gs
use_data(guideSetExampleFullAnnotation,
         compress="xz",
         overwrite=TRUE)




