context("Merge Peptides")
library(qPLEXanalyzer)

# data(human_anno)
# data(exp3_OHT_ESR1)
# MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#                                metadata = exp3_OHT_ESR1$metadata_qPLEX1,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")
MSnSet_phos <- readRDS("convertToMSnset_phospho_msnset.rds")
data(human_anno)

test_that("Merge peptides works", {
  expect_equal_to_reference(
    mergePeptides(MSnSet_data, sum, human_anno), 
    file="mergePeptides_msnset.rds")
  expect_equal_to_reference(
    mergePeptides(MSnSet_phos, sum, human_anno, 7), 
    file="mergePeptides_phos_msnset.rds")
})


