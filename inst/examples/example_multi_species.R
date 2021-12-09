library(crisprDesign)
library(crisprDesignDataS4)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10 
guides <- guideSetExample
guides <- addSpacerAlignments(guides,
                              txObject=txdb_human,
                              columnName="alnHuman")
guides <- guides[guides$n0==1]
guides <- addSpacerAlignments(guides,
                              txObject=txdb_mouse,
                              columnName="alnMouse",
                              bsgenome=bsgenome,
                              bowtie_index="mm10/mm10",
                              addSummary=FALSE)
onTargets(guides, "alnHuman")$cds
onTargets(guides, "alnMouse")$cds