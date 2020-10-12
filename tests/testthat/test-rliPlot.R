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

# test the check arg functions
test_that("argument checks - MSnset", {
    expect_error(rliPlot(MSnSetObj = 1),
                 regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - title", {
    expect_error(rliPlot(MSnSetObj = rawMSnSet, title = 1),
                 regexp = "title is not a string")
})
test_that("argument checks - sampleColours", {
    myCols <- assignColours(rawMSnSet)
    myCols[2] <- "blueq"
    expect_error(rliPlot(MSnSetObj = rawMSnSet,
                         sampleColours = myCols),
                 regexp = "sampleColours: blueq is not a valid colour")
    myCols <- c("red", "blue", "green", "pink")
    expect_error(rliPlot(MSnSetObj = rawMSnSet,
                         sampleColours = myCols),
                 regexp = "sampleColours should be a named vector")
    myCols <- assignColours(rawMSnSet)[-1]
    expect_error(rliPlot(MSnSetObj = rawMSnSet,
                         sampleColours = myCols),
                 regexp = "sampleColours should provide a colour for each")
    myCols <- assignColours(rawMSnSet)
    names(myCols)[1] <- "wibble"
    expect_error(rliPlot(MSnSetObj = rawMSnSet,
                         sampleColours = myCols),
                 regexp = "The names of sampleColours do not match the values")
})
test_that("argument checks - colourby", {
    expect_error(rliPlot(MSnSetObj = rawMSnSet, colourBy = 1),
                 regexp = "colourBy is not a string")
    expect_error(rliPlot(MSnSetObj = rawMSnSet, colourBy = "wibble"),
                 regexp = "column wibble not found in the MSnset metadata")
})
test_that("argument checks - omitIgG", {
    expect_error(rliPlot(MSnSetObj = rawMSnSet, omitIgG = 1),
                 regexp = "omitIgG is not a flag")
})
