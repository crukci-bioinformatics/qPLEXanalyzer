context("Assign Colours")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)

test_that("assign colours works", {
  expect_equal_to_reference(
    assignColours(rawMSnSet), 
    file="assignColours_sg.rds")
  expect_equal_to_reference(
    assignColours(rawMSnSet, colourBy = "BioRep"), 
    file="assignColours_rep.rds")
})
