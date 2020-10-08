context("Intensity density plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
exp3Met <- exp3_OHT_ESR1$metadata_qPLEX1[c(1,2,7,10),]
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[,c(1:6, c(1,2,7,10)+6)]
rawMSnSet <- convertToMSnset(exp3Int,
                               metadata = exp3Met,
                               indExpData = c(7:10),
                               Sequences = 2,
                               Accessions = 6)

# default settings
plt1 <- intensityPlot(rawMSnSet, title = "qPLEX_RIME_ER")

# colour by replicate
plt2 <- intensityPlot(rawMSnSet, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(rawMSnSet)$SampleGroup)))
names(customCols) <- unique(pData(rawMSnSet)$SampleGroup)
plt3 <- intensityPlot(rawMSnSet, title = "qPLEX_RIME_ER",
                        sampleColours = customCols)

test_that("intensityBoxplot works", {
  vdiffr::expect_doppelganger("Intensity density plot", plt1)
  vdiffr::expect_doppelganger("Intensity density plot colour by rep", plt2)
  vdiffr::expect_doppelganger("Intensity density plot custom colour", plt3)
})

# test argument checking

test_that("argument checks - MSnset", {
    expect_error(intensityPlot(MSnSetObj = 1),
                 regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - title", {
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, title = 1),
                 regexp = "title is not a string")
})
test_that("argument checks - sampleColours", {
    myCols <- assignColours(rawMSnSet)
    myCols[2] <- "blueq"
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "sampleColours: blueq is not a valid colour")
    myCols <- c("red", "blue", "green", "pink")
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "sampleColours should be a named vector")
    myCols <- assignColours(rawMSnSet)[-1]
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "sampleColours should provide a colour for each")
    myCols <- assignColours(rawMSnSet)
    names(myCols)[1] <- "wibble"
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "The names of sampleColours do not match the values")
})
test_that("argument checks - colourby", {
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  colourBy = 1),
                 regexp = "colourBy is not a string")
    expect_error(intensityPlot(MSnSetObj = rawMSnSet,
                                  colourBy = "wibble"),
                 regexp = "column wibble not found in the MSnset metadata")
})
test_that("argument checks - transform", {
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, transform = 1),
                 regexp = "transform is not a flag")
})
test_that("argument checks - xlab", {
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, xlab = 1),
                 regexp = "xlab is not a string")
})
test_that("argument checks - trFunc", {
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, trFunc = "mean"),
                 regexp = "trFunc is not a function")
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, trFunc = range),
                 regexp = "the specified function should tranform")
    expect_error(intensityPlot(MSnSetObj = rawMSnSet, trFunc = is.numeric),
                 regexp = "the specified function should tranform")
})
