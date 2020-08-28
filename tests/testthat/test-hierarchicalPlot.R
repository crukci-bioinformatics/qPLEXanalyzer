context("Hierarchical plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

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