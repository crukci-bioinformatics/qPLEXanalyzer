context("RLI plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:200, ]
rawMSnSet <- convertToMSnset(exp3Int,
                             metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)

# default
plt1 <- rliPlot(rawMSnSet, title = "qPLEX_RIME_ER")

# colour by rep and include IgG
plt2 <- rliPlot(rawMSnSet, omitIgG = FALSE, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(rawMSnSet)$SampleGroup)))
names(customCols) <- unique(pData(rawMSnSet)$SampleGroup)
plt3 <- rliPlot(rawMSnSet, sampleColours = customCols)

test_that("RLI plot works", {
  vdiffr::expect_doppelganger("RLI plot", plt1)
  vdiffr::expect_doppelganger("RLI plot colour by BioRep include IgG", plt2)
  vdiffr::expect_doppelganger("RLI plot custom colours", plt3)
})
