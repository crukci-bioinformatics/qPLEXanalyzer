context("Scaling normalizaion")
library(qPLEXanalyzer)

# data(exp2_Xlink)
# 
# MSnSet_data <- convertToMSnset(exp2_Xlink$intensities,
#                                metadata = exp2_Xlink$metadata,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)

MSnSet_data <- readRDS("convertToMSnset_exp2_msnset.rds")

test_that("Scaling normalization works", {
  expect_equal_to_reference(
    normalizeScaling(MSnSet_data, 
                     scalingFunction = median), 
    file="normalizeScaling_msnset.rds")
})

test_that("Scaling normalization with reference to protein works", {
  expect_equal_to_reference(
    normalizeScaling(MSnSet_data, 
                     scalingFunction = median, 
                     ProteinId = "P03372"), 
    file="normalizeScaling_protein_msnset.rds")
})