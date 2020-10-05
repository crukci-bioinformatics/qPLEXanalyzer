context("Compute differential analysis statistics")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1[1:1000,],
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)
MSnset_norm <- groupScaling(rawMSnSet, scalingFunction = median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", 
               tam.6h_vs_vehicle = "tam.6h - vehicle")
diffstats <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts)

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# # changes
# test_that("Compute Diff Stats works", {
#   expect_equal_to_reference(diffstats, file="computeDiffStats_msnset.rds")
# })

# The function creates:
# diffstats <- list(MSnSetObj = MSnSetObj, 
#                   fittedLM = fit, 
#                   fittedContrasts = fittedContrasts))
# Where the MSnSetObj is unchanged from input.
# Let's just compare the `fitted` objects


testObj <- list(fittedLM = diffstats$fittedLM,
                fittedContrasts = diffstats$fittedContrasts)

test_that("Compute Diff Stats works", {
  expect_equal_to_reference(testObj, file="computeDiffStats.rds")
})
