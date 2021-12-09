context('Testing findSpacers function')

library(crisprBase)
data(SpCas9, package="crisprBase")
data(enAsCas12a, package="crisprBase")
cas9 <- SpCas9
cas12a <- enAsCas12a

seq <- "CCAANAGTGAAACCACGTCTCTATAAAGAATACAAAAAATTAGCCGGGTGTTA"
results_findspacers_cas9 <- findSpacers(seq,
                                      crisprNuclease=cas9)
results_findspacers_cas9_notcanonical <- findSpacers(seq,
                                                   crisprNuclease=cas9,
                                                   canonical=FALSE)
results_findspacers_cas9_notcanonical_19 <- findSpacers(seq,
                                                      spacer_len=19,
                                                      crisprNuclease=cas9,
                                                      canonical=FALSE)
results_findspacers_cas12a <- findSpacers(seq,
                                        crisprNuclease=cas12a)
results_findspacers_cas12a_notcanonical <- findSpacers(seq,
                                                     crisprNuclease=cas12a,
                                                     canonical=FALSE)
results_findspacers_cas12a_notcanonical_22 <- findSpacers(seq,
                                                        spacer_len=22,
                                                        crisprNuclease=cas12a,
                                                        canonical=FALSE)

# test_that('Custom sequence', {
#     update=FALSE
#     expect_equal_to_reference(results_findspacers_cas9, update=update,
#         file=file.path("objects/results_findspacers_cas9.rds"))
#     expect_equal_to_reference(results_findspacers_cas9_notcanonical, update=update,
#         file=file.path("objects/results_findspacers_cas9_notcanonical.rds"))
#     expect_equal_to_reference(results_findspacers_cas9_notcanonical_19, update=update,
#         file=file.path("objects/results_findspacers_cas9_notcanonical_19.rds"))
#     expect_equal_to_reference(results_findspacers_cas12a, update=update,
#         file=file.path("objects/results_findspacers_cas12a.rds"))
#     expect_equal_to_reference(results_findspacers_cas12a_notcanonical, update=update,
#         file=file.path("objects/results_findspacers_cas12a_notcanonical.rds")) 
#     expect_equal_to_reference(results_findspacers_cas12a_notcanonical_22, update=update,
#         file=file.path("objects/results_findspacers_cas12a_notcanonical_22.rds"))  
# })




