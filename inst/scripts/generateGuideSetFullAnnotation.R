library(crisprDesign)
library(crisprDesignDataS4)
library(crisprScore)
library(biomaRt)
library(devtools)

txObject  <- txdb_human
tssObject <- tss_human
gr.repeats <- gr.repeats.hg38
bowtie_index <- "~/crisprIndices/bowtie/hg38/hg38"
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
                                   aligner='bowtie',
                                   txObject=txObject,
                                   tssObject=tssObject,
                                   n_mismatches=3,
                                   bowtie_index=bowtie_index)
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



