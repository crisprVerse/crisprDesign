library(crisprBase)
library(crisprDesignS4)
nuc <- CrisprNuclease(nucleaseName="JP",
                      metadata="My nuc",
                      pams="^N",
                      spacer_side="5prime",
                      spacer_length=5)
isCutting(nuc)
seq <- "AGGGGTTGGGTA"
gs <- findSpacers(seq,
                  both_strand=FALSE,
                  crisprNuclease=nuc)
gs <- addSpacerAlignments(gs,
                          n_mismatches=3,
                          custom_seq=seq,
                          aligner="biostrings")


### With a gap
nuc <- CrisprNuclease(nucleaseName="JP",
                      metadata="My nuc",
                      pams="^AGG",
                      spacer_gap=3,
                      spacer_side="5prime",
                      spacer_length=5)
hasSpacerGap(nuc)
seq <- "AGGGGTTGGGAGG"
gs <- findSpacers(seq,
                  both_strand=FALSE,
                  crisprNuclease=nuc)
gs <- addSpacerAlignments(gs,
                          n_mismatches=3,
                          custom_seq=seq,
                          aligner="biostrings")


