context("MA and Volcano plots")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1[1:1000,],
                             metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)
MSnset_norm <- groupScaling(rawMSnSet, scalingFunction = median)
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)

contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle", 
               tam.6h_vs_vehicle = "tam.6h - vehicle")
diffstats <- computeDiffStats(MSnSetObj=MSnset_Pnorm, contrasts=contrasts)

# MA plot
plt1 <- maVolPlot(diffstats, contrast = contrasts[1], plotType="MA")

# Volcano plot
plt2 <- maVolPlot(diffstats, contrast = contrasts[1], plotType="Volcano")

# Volcano plot with genes of interest "selected"
goi <- c("P48552", "Q96C36", "P54886", "P32322", "Q4ZG55", "P10644")
plt3 <- maVolPlot(diffstats, 
                  contrast = contrasts[1], 
                  plotType="Volcano", 
                  selectedGenes = goi)

test_that("MA plot works", {
  vdiffr::expect_doppelganger("MA plot", plt1)
})

test_that("Volcano plot works", {
  vdiffr::expect_doppelganger("Volcano plot", plt2)
})

test_that("Volcano plot with selected genes", {
    vdiffr::expect_doppelganger("Volcano plot with selected genes", plt3)
})
