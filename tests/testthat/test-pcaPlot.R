context("PCA plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

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
