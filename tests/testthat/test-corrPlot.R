context("Correlation plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

test_that("corrPlot works", {
  vdiffr::expect_doppelganger(
    "correlation plot default colour", 
    corrPlot(rawMSnSet, addValues=TRUE, title="Correlation plot"))
  vdiffr::expect_doppelganger(
    "correlation plot default colour no numbers", 
    corrPlot(rawMSnSet, addValues=FALSE, title="Correlation plot"))
  vdiffr::expect_doppelganger(
    "correlation plot yellow to pink", 
    corrPlot(rawMSnSet, addValues=TRUE, title="Correlation plot", 
             low_cor_colour="yellow", high_cor_colour="pink"))
  vdiffr::expect_doppelganger(
    "correlation plot change limit",
    corrPlot(rawMSnSet, addValues=TRUE, title="Correlation plot",
             low_cor_limit=0.5, high_cor_limit=1))
})
