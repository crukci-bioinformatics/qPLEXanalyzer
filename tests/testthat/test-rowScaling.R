context("Row scaling")
library(qPLEXanalyzer)

# data(exp3_OHT_ESR1)
# data(human_anno)
# 
# MSnset_reg <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX2,
#                               metadata = exp3_OHT_ESR1$metadata_qPLEX2,
#                               indExpData = c(7:16), 
#                               Sequences = 2, 
#                               Accessions = 6)
# MSnset_P <- summarizeIntensities(MSnset_reg, sum, human_anno)

MSnSet_P <- readRDS("summarizeIntensities_msnset.rds")

test_that("Row scaling works", {
  expect_equal_to_reference(
    rowScaling(MSnSet_P, mean), 
    file="rowScaling_msnset.rds")
})