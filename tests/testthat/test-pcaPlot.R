context("PCA plot")
library(qPLEXanalyzer)
data(exp3_OHT_ESR1)

rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

# default
exprs(rawMSnSet) <- exprs(rawMSnSet)+0.01
plt1 <- pcaPlot(rawMSnSet)

# colour by rep and omit IgG
plt2 <- pcaPlot(rawMSnSet, omitIgG = TRUE, colourBy = "BioRep")

# custom colours and PC2 v PC3
customCols <- rainbow(length(unique(pData(rawMSnSet)$SampleGroup)))
names(customCols) <- unique(pData(rawMSnSet)$SampleGroup)
plt3 <- pcaPlot(rawMSnSet, sampleColours = customCols, x.PC=2)

test_that("PCA plot works", {
  vdiffr::expect_doppelganger("PCA plot", plt1)
  vdiffr::expect_doppelganger("PCA plot colour by BioRep omit IgG", plt2)
  vdiffr::expect_doppelganger("PCA plot custom colours and PC2 v PC3", plt3)
})

