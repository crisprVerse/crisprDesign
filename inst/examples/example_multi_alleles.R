library(crisprDesign)
library(crisprDesignDataS4)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
bsgenome_major <- BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major
bsgenome_minor <- BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor


guides <- guideSetExample
n_mismatches=0
guides <- addSpacerAlignments(guides,
                              n_mismatches=n_mismatches,
                              txObject=txdb_human,
                              bsgenome=bsgenome_major,
                              bowtie_index="hg38_major/hg38_major",
                              columnName="alnMajor")
guides <- guides[guides$n0==1]
guides <- addSpacerAlignments(guides,
                              n_mismatches=n_mismatches,
                              addSummary=FALSE,
                              txObject=txdb_human,
                              bsgenome=bsgenome_minor,
                              bowtie_index="hg38_minor/hg38_minor",
                              columnName="alnMinor")

a=onTargets(guides, "alnMajor")
b=onTargets(guides, "alnMinor")
length(a)
length(b)





