context("Quantile normalizaion")
library(qPLEXanalyzer)

# data(exp2_Xlink)
# 
# MSnSet_data <- convertToMSnset(exp2_Xlink$intensities,
#                                metadata = exp2_Xlink$metadata,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)

MSnSet_data <- readRDS("convertToMSnset_exp2_msnset.rds")

test_that("Quantile normalization works", {
  expect_equal_to_reference(
    normalizeQuantiles(MSnSet_data), 
    file="normalizeQuantiles_msnset.rds")
})