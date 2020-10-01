context("MA and Volcano plots")
library(qPLEXanalyzer)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
diffstats <- readRDS("computeDiffStats_msnset.rds")

plt1 <- maVolPlot(diffstats, contrast = contrasts, plotType="MA")
plt2 <- maVolPlot(diffstats, contrast = contrasts, plotType="Volcano")
goi <- c("P48552", "Q96C36", "P54886", "P32322", "Q4ZG55", "P10644")
plt3 <- maVolPlot(diffstats, 
                  contrast = contrasts, 
                  plotType="Volcano", 
                  selectedGenes = goi)

test_that("MA plot works", {
  vdiffr::expect_doppelganger("MA plot", plt1)
})


test_that("Volcano plot works", {
  vdiffr::expect_doppelganger("Volcano plot", plt2)
})


test_that("Volcano plot with selected genes", {
    vdiffr::expect_doppelganger("Volcano plot with selected genes", plt3)
})
