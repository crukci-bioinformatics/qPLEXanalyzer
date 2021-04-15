context("Group Scaling normalization")
library(qPLEXanalyzer)

data(exp2_Xlink)

exp2Int <- exp2_Xlink$intensities[1:500, ]
rawMSnSet <- convertToMSnset(exp2Int,
                             metadata = exp2_Xlink$metadata,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
normMSnSet <- groupScaling(rawMSnSet)

# Test the function

# The function only changes the expression set, so let's just compare that

testObj <-  exprs(normMSnSet)

test_that("Group scaling normalization works", {
    expect_equal_to_reference(testObj, 
                              file="groupScaling.rds")
})

# Test the argument checks

test_that("argument checks - MSnset", {
    expect_error(groupScaling(1),
                 regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - scalingFunction", {
    expect_error(groupScaling(rawMSnSet, scalingFunction = sum), 
                 regexp = "scalingFunction should be mean or median")
})
test_that("argument checks - groupingColumn", {
    expect_error(groupScaling(rawMSnSet, groupingColumn = 1), 
                 regexp = "groupingColumn is not a string")
    expect_error(groupScaling(rawMSnSet, groupingColumn = "Wibble"), 
                 regexp = "column Wibble not found in the MSnset metadata")
})
