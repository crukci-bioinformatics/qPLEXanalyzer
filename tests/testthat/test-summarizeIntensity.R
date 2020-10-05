context("Summarize intensities to protein level")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:2000, ]
rawMSnSet <- convertToMSnset(exp3Int,
                             metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)

protMSnSet <- summarizeIntensities(rawMSnSet, sum, human_anno)


# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# # changes
# test_that("Summarize intensities works", {
#   expect_equal_to_reference(protMSnSet, 
#                         file="summarizeIntensities_msnset.rds")
# })

# The function creates an entirely new MSnSet obj, we should compare the samples
# features data, and merged intensities of each object

protTestList <- list(Samples = pData(protMSnSet),
                     Features = fData(protMSnSet),
                     MergedIntensitied = exprs(protMSnSet))

test_that("Summarize intensities works", {
    expect_equal_to_reference(protTestList, 
                              file="summarizeIntensities.rds")
})