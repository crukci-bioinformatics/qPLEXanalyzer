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

# test check argument errors

test_that("argument checks - diffstats", {
    expect_error(maVolPlot(diffstats = 1, contrast = contrasts[1]),
                 regexp = "diffstats is not a valid output")
})
test_that("argument checks - contrast", {
    expect_error(maVolPlot(diffstats = diffstats, contrast = 1),
                 regexp = "contrast is not a string")
    expect_error(maVolPlot(diffstats = diffstats, contrast = "wibble"),
                 regexp = "'wibble' is not a valid contrast")
})
test_that("argument checks - title", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           title = 1),
                 regexp = "title is not a string")
})
test_that("argument checks - controlGroup", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           controlGroup = 1),
                 regexp = "controlGroup has to be a string .* or NULL")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           controlGroup = "Wibble"),
                 regexp = "controlGroup: 'Wibble' is not found in .* diffstats")
})
test_that("argument checks - selectedGenes", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           selectedGenes = 1),
                 regexp = "'selectedGenes' should be .* of Accessions or NULL")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           selectedGenes = c("A0AV96", "Wibble")),
                 regexp = "Some of the genes provided in 'selectedGenes'")
})
test_that("argument checks - fdrCutOff", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           fdrCutOff = "A"),
                 regexp = "fdrCutOff is not a number")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           fdrCutOff = c(0.2, 0.5)),
                 regexp = "fdrCutOff .*(a length one numeric vector)")
})
test_that("argument checks - lfcCutOff", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           lfcCutOff = "A"),
                 regexp = "lfcCutOff is not a number")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           lfcCutOff = c(0.2, 0.5)),
                 regexp = "lfcCutOff .*(a length one numeric vector)")
})
test_that("argument checks - controlLfcCutOff", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           controlLfcCutOff = "A"),
                 regexp = "controlLfcCutOff is not a number")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           controlLfcCutOff = c(0.2, 0.5)),
                 regexp = "controlLfcCutOff .*(a length one numeric vector)")
})
test_that("argument checks - plotType", {
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           plotType = 1),
                 regexp = "plotType is not a string")
    expect_error(maVolPlot(diffstats = diffstats,
                           contrast = contrasts[1],
                           plotType = c("MA", "Vol")),
                 regexp = "plotType .*(a length one character vector)")
})
