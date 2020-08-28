context("Intensity density plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

# default settings
plt1 <- intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# colour by replicate
plt2 <- intensityBoxplot(MSnSet_data, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
plt3 <- intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER",
                        sampleColours = customCols)

test_that("intensityBoxplot works", {
  vdiffr::expect_doppelganger("Intensity density plot", plt1)
  vdiffr::expect_doppelganger("Intensity density plot colour by rep", plt2)
  vdiffr::expect_doppelganger("Intensity density plot custom colour", plt3)
})