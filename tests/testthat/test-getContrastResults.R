context("Get contrast results table")
library(qPLEXanalyzer)

# data(human_anno)
# data(exp3_OHT_ESR1)
# 
# MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#                                metadata = exp3_OHT_ESR1$metadata_qPLEX1,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)
# MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
# MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
# diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")

diffstats <- readRDS("computeDiffStats_msnset.rds")

test_that("Get Contrast results works", {
  expect_equal_to_reference(
    getContrastResults(diffstats=diffstats, contrast=contrasts), 
    file="getContrastResults.rds")
})