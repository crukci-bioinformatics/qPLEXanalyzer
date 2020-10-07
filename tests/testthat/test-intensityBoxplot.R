context("Intensity boxplot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:200, ]
rawMSnSet <- convertToMSnset(exp3Int,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

# default settings
plt1 <- intensityBoxplot(rawMSnSet, title = "qPLEX_RIME_ER")

# colour by replicate
plt2 <- intensityBoxplot(rawMSnSet, colourBy = "BioRep")

# custom colours
customCols <- rainbow(length(unique(pData(rawMSnSet)$SampleGroup)))
names(customCols) <- unique(pData(rawMSnSet)$SampleGroup)
plt3 <- intensityBoxplot(rawMSnSet, 
                         sampleColours = customCols)

test_that("intensityBoxplot works", {
  vdiffr::expect_doppelganger("Intensity boxplot", plt1)
  vdiffr::expect_doppelganger("Intensity boxplot colour by rep", plt2)
  vdiffr::expect_doppelganger("Intensity boxplot custom colour", plt3)
})


# the argument checks
test_that("argument checks - MSnset", {
    expect_error(intensityBoxplot(MSnSetObj = 1), 
                 regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - title", {
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet, title = 1), 
                 regexp = "title is not a string")
})
test_that("argument checks - sampleColours", {
    myCols <- assignColours(rawMSnSet)
    myCols[2] <- "blueq"
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols), 
                 regexp = "sampleColours: blueq is not a valid colour")
    myCols <- c("red", "blue", "green", "pink")
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),  
                 regexp = "sampleColours should be a named vector")
    myCols <- assignColours(rawMSnSet)[-1]
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "sampleColours should provide a colour for each")
    myCols <- assignColours(rawMSnSet)
    names(myCols)[1] <- "wibble"
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet,
                                  sampleColours = myCols),
                 regexp = "The names of sampleColours do not match the values")
})
test_that("argument checks - colourby", {
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet, 
                                  colourBy = 1), 
                 regexp = "colourBy is not a string")
    expect_error(intensityBoxplot(MSnSetObj = rawMSnSet, 
                                  colourBy = "wibble"), 
                 regexp = "column wibble not found in the MSnset metadata")
})




