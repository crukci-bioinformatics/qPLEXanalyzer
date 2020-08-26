context("Regression analysis")
library(qPLEXanalyzer)

# data(exp3_OHT_ESR1)
# data(human_anno)
# 
# MSnset_reg <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX2,
#                               metadata = exp3_OHT_ESR1$metadata_qPLEX2,
#                               indExpData = c(7:16), 
#                               Sequences = 2, 
#                               Accessions = 6)
# MSnset_P <- summarizeIntensities(MSnset_reg, sum, human_anno)
# MSnset_P <- rowScaling(MSnset_P, mean)

MSnset_P <- readRDS("rowScaling_msnset.rds")
IgG_ind <- which(pData(MSnset_P)$SampleGroup == "IgG")

test_that("Regeression analysis works", {
  expect_equal_to_reference(
    regressIntensity(MSnset_P, controlInd = IgG_ind, ProteinId = "P03372"), 
    file="regressIntensities_msnset.rds")
})