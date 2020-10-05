context("Quantile normalizaion")
library(qPLEXanalyzer)
data(exp2_Xlink)

exp2Int <- exp2_Xlink$intensities[1:1000, ]
rawMSnSet <- convertToMSnset(exp2Int,
                             metadata = exp2_Xlink$metadata,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
normMSnSet <- normalizeQuantiles(rawMSnSet)

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# test_that("Quantile normalization works", {
#   expect_equal_to_reference(normMSnSet, file="normalizeQuantiles_msnset.rds")
# })

# The function only changes the expression set, so let's just compare that

testObj <-  exprs(normMSnSet)

test_that("Quantile normalization works", {
    expect_equal_to_reference(testObj, 
                              file="normalizeQuantiles.rds")
})
