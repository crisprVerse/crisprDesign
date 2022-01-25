library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
tx_id <- "ENST00000538872"
df <- getTxInfoDataFrame(tx_id=tx_id,
                         txObject=ensembl_human,
                         bsgenome=bsgenome)
df <- addEditedAlleles(df)
plot(df$)
