context("Correlation plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

test_that("corrPlot works", {
  vdiffr::expect_doppelganger(
    "correlation plot default colour", 
    corrPlot(MSnSet_data, addValues=TRUE, title="Correlation plot"))
  vdiffr::expect_doppelganger(
    "correlation plot default colour no numbers", 
    corrPlot(MSnSet_data, addValues=FALSE, title="Correlation plot"))
  vdiffr::expect_doppelganger(
    "correlation plot yellow to pink", 
    corrPlot(MSnSet_data, addValues=TRUE, title="Correlation plot", 
             low_cor_colour="yellow", high_cor_colour="pink"))
})