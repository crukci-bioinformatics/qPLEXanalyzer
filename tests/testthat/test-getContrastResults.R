context("Get contrast results table")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1[1:2000,],
                             metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
#MSnset_norm <- groupScaling(rawMSnSet, scalingFunction = median)
#MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
MSnset_Pnorm <- summarizeIntensities(rawMSnSet, sum, human_anno)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", 
               tam.6h_vs_vehicle = "tam.6h - vehicle")
diffstats <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")

res <- getContrastResults(diffstats=diffstats, contrast=contrasts)

# test the function

test_that("Get Contrast results works", {
  expect_equal_to_reference(res, file="getContrastResults.rds")
})

# test the argument checks
test_that("argument checks - diffstats", {
    testDat <- list(MSnSetObj = 1, fittedLM = 1, Wibble = 1)
    errMsg <- "diffstats is not a valid output of the computeDiffStats function"
    expect_error(getContrastResults(diffstats=testDat, contrast=contrasts), 
                 regexp = errMsg)
})
test_that("argument checks - contrast", {
    expect_error(getContrastResults(diffstats=diffstats, contrast=1), 
                 regexp = "contrast is not a string")
    expect_error(getContrastResults(diffstats=diffstats, contrast = "Wibble"), 
               regexp = "'Wibble' is not a valid contrast. Available contrasts")
})
test_that("argument checks - controlGroup", {
    expect_error(getContrastResults(diffstats=diffstats, 
                                    contrast=contrasts,
                                    controlGroup = 1), 
                 regexp = "controlGroup has to be a string")
    expect_error(getContrastResults(diffstats=diffstats, 
                                    contrast = contrasts,
                                    controlGroup = "Wibble"),  
         regexp = "controlGroup: 'Wibble' is not found in the diffstats object")
})