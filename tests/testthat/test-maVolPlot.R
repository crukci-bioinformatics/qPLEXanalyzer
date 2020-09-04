context("MA and Volcano plots")
library(qPLEXanalyzer)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- readRDS("computeDiffStats_msnset.rds")

plt1 <- maVolPlot(diffstats, contrast = contrasts, plotType="MA")
plt2 <- maVolPlot(diffstats, contrast = contrasts, plotType="Volcano")

test_that("MA plot works", {
  vdiffr::expect_doppelganger("MA plot", plt1)
})


test_that("Volcano plot works", {
  vdiffr::expect_doppelganger("Volcano plot", plt2)
})
