context("RLI plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

rliPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
rliPlot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)

# default
plt1 <- rliPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# colour by rep and include IgG
plt2 <- rliPlot(MSnSet_data, omitIgG = FALSE, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
plt3 <- rliPlot(MSnSet_data, sampleColours = customCols)

test_that("RLI plot works", {
  vdiffr::expect_doppelganger("RLI plot", plt1)
  vdiffr::expect_doppelganger("RLI plot colour by BioRep include IgG", plt2)
  vdiffr::expect_doppelganger("RLI plot custom colours", plt3)
})
