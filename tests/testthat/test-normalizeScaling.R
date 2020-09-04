context("Scaling normalizaion")
library(qPLEXanalyzer)
data(exp2_Xlink)

exp2Int <- exp2_Xlink$intensities[1:50, ]
rawMSnSet <- convertToMSnset(exp2Int,
                             metadata = exp2_Xlink$metadata,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
normMSnSet <- normalizeScaling(rawMSnSet, scalingFunction = median)
pnormMSnSet <- normalizeScaling(rawMSnSet, 
                                scalingFunction = median, 
                                ProteinId = "P04264")

test_that("Scaling normalization works", {
  expect_equal_to_reference(normMSnSet, file="normalizeScaling_msnset.rds")
})

test_that("Scaling normalization with reference to protein works", {
  expect_equal_to_reference(pnormMSnSet, file="normalizeScaling_prot_msnset.rds")
})