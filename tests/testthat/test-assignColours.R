context("Assign Colours")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")

test_that("assign colours works", {
  expect_equal_to_reference(
    assignColours(MSnSet_data), 
    file="assignColours_sg.rds")
  expect_equal_to_reference(
    assignColours(MSnSet_data, colourBy = "BioRep"), 
    file="assignColours_rep.rds")
})