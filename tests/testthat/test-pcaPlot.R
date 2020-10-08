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

# test the check arg functions
test_that("argument checks - MSnset", {
    expect_error(pcaPlot(MSnSetObj = 1),
                 regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - omitIgG", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, omitIgG = 1),
               regexp = "omitIgG is not a flag")
})
test_that("argument checks - omitIgG", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, omitIgG = 1),
               regexp = "omitIgG is not a flag")
})
test_that("argument checks - sampleColours", {
  myCols <- assignColours(rawMSnSet)
  myCols[2] <- "blueq"
  expect_error(pcaPlot(MSnSetObj = rawMSnSet,
                       sampleColours = myCols),
               regexp = "sampleColours: blueq is not a valid colour")
  myCols <- c("red", "blue", "green", "pink")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet,
                       sampleColours = myCols),
               regexp = "sampleColours should be a named vector")
  myCols <- assignColours(rawMSnSet)[-1]
  expect_error(pcaPlot(MSnSetObj = rawMSnSet,
                       sampleColours = myCols),
               regexp = "sampleColours should provide a colour for each")
  myCols <- assignColours(rawMSnSet)
  names(myCols)[1] <- "wibble"
  expect_error(pcaPlot(MSnSetObj = rawMSnSet,
                       sampleColours = myCols),
               regexp = "The names of sampleColours do not match the values")
})
test_that("argument checks - transFunc", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, transFunc = "mean"),
               regexp = "transFunc is not a function")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, transFunc = range),
               regexp = "the specified function should tranform")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, transFunc = is.numeric),
               regexp = "the specified function should tranform")
})
test_that("argument checks - transform", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, transform = 1),
               regexp = "transform is not a flag")
})
test_that("argument checks - colourby", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, colourBy = 1),
               regexp = "colourBy is not a string")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, colourBy = "wibble"),
               regexp = "column wibble not found in the MSnset metadata")
})
test_that("argument checks - title", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, title = 1),
               regexp = "title is not a string")
})
test_that("argument checks - labelColumn", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, labelColumn = 1),
               regexp = "labelColumn is not a string")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, labelColumn = "Wibble"),
               regexp = "column Wibble not found in the MSnset metadata")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, 
                       labelColumn = c("SampleGroup", "BioRep")),
               regexp = "labelColumn is not a string")
})
test_that("argument checks - labelsize", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, labelsize = "A"),
               regexp = "labelsize is not a number")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, labelsize = c(1, 1.5)),
               regexp = "labelsize is not a number")
})
test_that("argument checks - pointsize", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, pointsize = "A"),
               regexp = "pointsize is not a number")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, pointsize = c(1, 1.5)),
               regexp = "pointsize is not a number")
})
test_that("argument checks - x.nudge", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, x.nudge = "A"),
               regexp = "x.nudge is not a number")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, x.nudge = c(1, 1.5)),
               regexp = "x.nudge is not a number")
})
test_that("argument checks - x.PC", {
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, x.PC = "A"),
               regexp = "x.PC is not a count")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, x.PC = 1.5),
               regexp = "x.PC is not a count")
  expect_error(pcaPlot(MSnSetObj = rawMSnSet, x.PC = c(1, 1.5)),
               regexp = "x.PC is not a count")
})
