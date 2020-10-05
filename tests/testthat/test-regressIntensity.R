context("Regression analysis")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
data(human_anno)

rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX2[1:2000,],
                              metadata = exp3_OHT_ESR1$metadata_qPLEX2,
                              indExpData = c(7:16),
                              Sequences = 2,
                              Accessions = 6)
MSnset_P <- summarizeIntensities(rawMSnSet, sum, human_anno)
MSnset_P <- rowScaling(MSnset_P, mean)

IgG_ind <- which(pData(MSnset_P)$SampleGroup == "IgG")

regrMSnSet <- regressIntensity(MSnset_P, 
                                  controlInd = IgG_ind, 
                                  ProteinId = "A0AV96", 
                                  plot=FALSE)

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# # changes
# test_that("Regeression analysis works", {
#   expect_equal_to_reference(
#     regressIntensity(MSnset_P, controlInd = IgG_ind, ProteinId = "A0AV96", plot=FALSE), 
#     file="regressIntensities_msnset.rds")
# })

# the function updates the expresssion set and the sample data, let's test these

testList <- list(Samples = pData(regrMSnSet),
                 RegInt = exprs(regrMSnSet))

test_that("Regeression analysis works", {
    expect_equal_to_reference(testList, file="regressIntensities.rds")
})
