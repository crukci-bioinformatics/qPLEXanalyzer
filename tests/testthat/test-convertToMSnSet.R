context("Convert to MSnSet")
library(qPLEXanalyzer)
data("exp2_Xlink")

test_that("Convert works", {
  expect_equal_to_reference(
    convertToMSnset(exp2_Xlink$intensities,
                  metadata = exp2_Xlink$metadata,
                  indExpData = c(7:16), 
                  Sequences = 2, 
                  Accessions = 6), 
    file="convertToMSnset_exp2_msnset.rds")
  
  expect_equal_to_reference(
    convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                    metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                    indExpData = c(7:16), 
                    Sequences = 2, 
                    Accessions = 6), 
    file="convertToMSnset_oht_esr1_msnset.rds")
})