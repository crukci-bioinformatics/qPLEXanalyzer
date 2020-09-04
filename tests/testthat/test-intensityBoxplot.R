context("Intensity boxplot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:200, ]
MSnSet_data <- convertToMSnset(exp3Int,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

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
