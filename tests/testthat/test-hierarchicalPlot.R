context("Hierarchical plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

# Default
plt1 <- hierarchicalPlot(MSnSet_data, title = "qPLEX_RIME_ER")
# Colour by rep
plt2 <- hierarchicalPlot(MSnSet_data, colourBy = "BioRep")
# Custom colours
myCols <- c("green", "blue", "red", "black", "pink") %>% 
  setNames(unique(pData(MSnSet_data)$SampleGroup))
plt3 <- hierarchicalPlot(MSnSet_data, sampleColours=myCols)
# Vertical version
plt4 <- hierarchicalPlot(MSnSet_data, horizontal = FALSE)

test_that("hierarchicalPlot works", {
  vdiffr::expect_doppelganger("Hierarchical Plot", plt1)
  vdiffr::expect_doppelganger("Hierarchical Plot colour by rep", plt2)
  vdiffr::expect_doppelganger("Hierarchical Plot custom colour", plt3)
  vdiffr::expect_doppelganger("Hierarchical Plot veritcal", plt4)
})
