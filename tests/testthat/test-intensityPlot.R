context("Intensity density plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
exp3Met <- exp3_OHT_ESR1$metadata_qPLEX1[c(1,2,7,10),]
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[,c(1:6, c(1,2,7,10)+6)]
MSnSet_data <- convertToMSnset(exp3Int,
                               metadata = exp3Met,
                               indExpData = c(7:10),
                               Sequences = 2,
                               Accessions = 6)

# default settings
plt1 <- intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# colour by replicate
plt2 <- intensityPlot(MSnSet_data, colourBy = "BioRep")

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