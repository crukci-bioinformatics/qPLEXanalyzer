context("PCA plot")
library(qPLEXanalyzer)
data(exp3_OHT_ESR1)

MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

# default
exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
plt1 <- pcaPlot(MSnSet_data)

# colour by rep and omit IgG
plt2 <- pcaPlot(MSnSet_data, omitIgG = TRUE, colourBy = "BioRep")

# custom colours and PC2 v PC3
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
plt3 <- pcaPlot(MSnSet_data, sampleColours = customCols, x.PC=2)

test_that("PCA plot works", {
  vdiffr::expect_doppelganger("PCA plot", plt1)
  vdiffr::expect_doppelganger("PCA plot colour by BioRep omit IgG", plt2)
  vdiffr::expect_doppelganger("PCA plot custom colours and PC2 v PC3", plt3)
})

