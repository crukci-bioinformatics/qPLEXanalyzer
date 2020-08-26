context("Group Scaling normalizaion")
library(qPLEXanalyzer)

# data(exp2_Xlink)
# 
# MSnSet_data <- convertToMSnset(exp2_Xlink$intensities,
#                                metadata = exp2_Xlink$metadata,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)

MSnSet_data <- readRDS("convertToMSnset_exp2_msnset.rds")

test_that("Group scaling normalization works", {
  expect_equal_to_reference(
    groupScaling(MSnSet_data), 
    file="groupScaling_msnset.rds")
})