context("Get contrast results table")
library(qPLEXanalyzer)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")

diffstats <- readRDS("computeDiffStats_msnset.rds")

test_that("Get Contrast results works", {
  expect_equal_to_reference(
    getContrastResults(diffstats=diffstats, contrast=contrasts), 
    file="getContrastResults.rds")
})
