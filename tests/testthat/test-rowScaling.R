context("Row scaling")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
data(human_anno)

MSnset_reg <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX2[1:1000,],
                              metadata = exp3_OHT_ESR1$metadata_qPLEX2,
                              indExpData = c(7:16),
                              Sequences = 2,
                              Accessions = 6)
MSnSet_P <- summarizeIntensities(MSnset_reg, sum, human_anno)

rnormMSnSet <- rowScaling(MSnSet_P, mean)

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# # changes
# test_that("Row scaling works", {
#     expect_equal_to_reference(
#         rowScaling(MSnSet_P, mean), 
#         file="rowScaling_msnset.rds")
# })

# The function only changes the expression set, so let's just compare that

testObj <-  exprs(rnormMSnSet)

test_that("Row scaling works", {
    expect_equal_to_reference(testObj, file="rowScaling.rds")
})
