context("Scaling normalizaion")
library(qPLEXanalyzer)
data(exp2_Xlink)

exp2Int <- exp2_Xlink$intensities[1:1000, ]
rawMSnSet <- convertToMSnset(exp2Int,
                             metadata = exp2_Xlink$metadata,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
normMSnSet <- normalizeScaling(rawMSnSet, scalingFunction = median)
pnormMSnSet <- normalizeScaling(rawMSnSet, 
                                scalingFunction = median, 
                                ProteinId = "P04264")

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# test_that("Scaling normalization works", {
#   expect_equal_to_reference(normMSnSet, file="normalizeScaling_msnset.rds")
# })
# 
# test_that("Scaling normalization with reference to protein works", {
#   expect_equal_to_reference(pnormMSnSet, file="normalizeScaling_prot_msnset.rds")
# })


# The function only changes the expression set, so let's just compare that

testObj <-  exprs(normMSnSet)

test_that("Scaling normalization works", {
    expect_equal_to_reference(testObj, file="normalizeScaling.rds")
})

ptestObj <-  exprs(pnormMSnSet)

test_that("Scaling normalization with reference to protein works", {
    expect_equal_to_reference(ptestObj, file="normalizeScaling_prot.rds")
})
