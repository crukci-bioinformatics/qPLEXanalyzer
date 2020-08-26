context("Compute differential analysis statistics")
library(qPLEXanalyzer)

# data(human_anno)
# data(exp3_OHT_ESR1)
# MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
#                                metadata = exp3_OHT_ESR1$metadata_qPLEX1,
#                                indExpData = c(7:16), 
#                                Sequences = 2, 
#                                Accessions = 6)
# MSnset_norm <- groupScaling(MSnSet_data, scalingFunction = median)
# MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)

MSnset_Pnorm <- readRDS("summarizeIntensities_msnset.rds")
contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", 
               tam.6h_vs_vehicle = "tam.6h - vehicle")


test_that("Compute Diff Stats works", {
  expect_equal_to_reference(
    computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts), 
    file="computeDiffStats_msnset.rds")
})