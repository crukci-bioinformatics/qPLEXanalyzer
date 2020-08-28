context("Intensity boxplot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

# default settings
plt1 <- intensityBoxplot(MSnSet_data, title = "qPLEX_RIME_ER")

# colour by replicate
plt2 <- intensityBoxplot(MSnSet_data, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
plt3 <- intensityBoxplot(MSnSet_data, 
                         sampleColours = customCols)

test_that("intensityBoxplot works", {
  vdiffr::expect_doppelganger("Intensity boxplot", plt1)
  vdiffr::expect_doppelganger("Intensity boxplot colour by rep", plt2)
  vdiffr::expect_doppelganger("Intensity boxplot custom colour", plt3)
})