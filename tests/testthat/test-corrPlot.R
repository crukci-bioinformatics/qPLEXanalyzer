context("Correlation plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

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
