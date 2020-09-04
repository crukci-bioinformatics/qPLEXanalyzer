context("Quantile normalizaion")
library(qPLEXanalyzer)
data(exp2_Xlink)

exp2Int <- exp2_Xlink$intensities[1:50, ]
rawMSnSet <- convertToMSnset(exp2Int,
                             metadata = exp2_Xlink$metadata,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
normMSnSet <- normalizeQuantiles(rawMSnSet)


test_that("Quantile normalization works", {
  expect_equal_to_reference(normMSnSet, file="normalizeQuantiles_msnset.rds")
})
