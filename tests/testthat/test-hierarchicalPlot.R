context("Hierarchical plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

# Default
plt1 <- hierarchicalPlot(rawMSnSet, title = "qPLEX_RIME_ER")
# Colour by rep
plt2 <- hierarchicalPlot(rawMSnSet, colourBy = "BioRep")
# Custom colours
myCols <- c("green", "blue", "red", "black", "pink") %>% 
  setNames(unique(pData(rawMSnSet)$SampleGroup))
plt3 <- hierarchicalPlot(rawMSnSet, sampleColours=myCols)
# Vertical version
plt4 <- hierarchicalPlot(rawMSnSet, horizontal = FALSE)

# test the function

test_that("hierarchicalPlot works", {
  vdiffr::expect_doppelganger("Hierarchical Plot", plt1)
  vdiffr::expect_doppelganger("Hierarchical Plot colour by rep", plt2)
  vdiffr::expect_doppelganger("Hierarchical Plot custom colour", plt3)
  vdiffr::expect_doppelganger("Hierarchical Plot veritcal", plt4)
})

# test the argument checks

test_that("argument checks - MSnset", {
  expect_error(hierarchicalPlot(MSnSetObj = 1),
               regexp = "MSnSetObj has to be of class MSnSet")
})
test_that("argument checks - colourBy", {
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, colourBy = 1), 
               regexp = "colourBy is not a string")
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, colourBy = "Wibble"), 
               regexp = "column Wibble not found in the MSnset metadata")
})
test_that("argument checks - sampleColours", {
  myCols <- assignColours(rawMSnSet)
  myCols[2] <- "blueq"
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, sampleColours = myCols),
               regexp = "sampleColours: blueq is not a valid colour")
  myCols <- c("red", "blue", "green", "pink")
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, sampleColours = myCols),
               regexp = "sampleColours should be a named vector")
  myCols <- assignColours(rawMSnSet)[-1]
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, sampleColours = myCols),
               regexp = "sampleColours should provide a colour for each")
  myCols <- assignColours(rawMSnSet)
  names(myCols)[1] <- "wibble"
  expect_error(hierarchicalPlot(MSnSetObj = rawMSnSet, sampleColours = myCols),
               regexp = "The names of sampleColours do not match the values")
})
