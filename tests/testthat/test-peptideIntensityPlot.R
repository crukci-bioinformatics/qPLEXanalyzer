context("Peptide intensity plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
data(human_anno)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1
exp3Int <- exp3Int[exp3Int$Master.Protein.Accessions=="P03372",]
rawMSnSet <-convertToMSnset(exp3Int,
                            metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                            indExpData = c(7:16),
                            Sequences = 2,
                            Accessions = 6)
protMSnSet <- summarizeIntensities(rawMSnSet, sum, human_anno)

# standard
plt1 <- peptideIntensityPlot(rawMSnSet, 
                     combinedIntensities=protMSnSet, 
                     ProteinID="P03372", 
                     ProteinName= "ESR1")

# no combined intensity
plt2 <- peptideIntensityPlot(rawMSnSet, 
                           ProteinID="P03372", 
                           ProteinName= "ESR1")

# selected sequence
plt3 <- peptideIntensityPlot(rawMSnSet, 
                             combinedIntensities=protMSnSet,
                             ProteinID="P03372", 
                             ProteinName= "ESR1",
                             selectedSequence = "[K].NVVPLYDLLLEMLDAHR.[L]")

# selected sequence & modification
plt4 <- peptideIntensityPlot(rawMSnSet, 
                             combinedIntensities=protMSnSet,
                             ProteinID="P03372", 
                             ProteinName= "ESR1",
                             selectedSequence = "[K].NVVPLYDLLLEMLDAHR.[L]",
                             selectedModifications = "1xTMT6plex [N-Term]")


test_that("Peptide plot works", {
  vdiffr::expect_doppelganger("Standard", plt1)
  vdiffr::expect_doppelganger("No combined intensity", plt2)
  vdiffr::expect_doppelganger("With selected sequence", plt3)
  vdiffr::expect_doppelganger("With selected modification", plt4)
})

# test argument checks

test_that("argument checks - MSnset", {
    expect_error(peptideIntensityPlot(MSnSetObj = 1, 
                                      ProteinID = "P03372", 
                                      ProteinName = "ESR1"),
                 regexp = "MSnSetObj has to be of class MSnSet")
    expect_error(peptideIntensityPlot(MSnSetObj = protMSnSet,
                                      ProteinID= "P03372", 
                                      ProteinName = "ESR1"),
                 regexp = "This MSnSet is not a peptide data set")
})
test_that("argument checks - MSnset", {
    expect_error(peptideIntensityPlot(MSnSetObj = 1,
                                      ProteinID = "P03372", 
                                      ProteinName = "ESR1"),
                 regexp = "MSnSetObj has to be of class MSnSet")
    expect_error(peptideIntensityPlot(MSnSetObj = protMSnSet, 
                                      ProteinID = "P03372", 
                                      ProteinName = "ESR1"),
                 regexp = "This MSnSet is not a peptide data set")
})
test_that("argument checks - ProteinId", {
    expect_error(peptideIntensityPlot(MSnSetObj = rawMSnSet, 
                                      ProteinID = 1, 
                                      ProteinName = "ESR1"),
                 regexp = "ProteinID is not a string")
    expect_error(peptideIntensityPlot(MSnSetObj = rawMSnSet, 
                                      ProteinID = c("P03372", "A889hi9"), 
                                      ProteinName = "ESR1"),
                 regexp = "ProteinID is not a string")
    expect_error(peptideIntensityPlot(MSnSetObj = rawMSnSet, 
                                      ProteinID = "P03372xxxx", 
                                      ProteinName = "ESR1"),
                 regexp = "ProteinID .* is not found the MSnset feature data")
})
test_that("argument checks - combinedIntensities", {
    errMsg <- str_c("combinedIntensities should be either NULL or an object ",
                    "of class MSnSet containing summarized protein data")
    expect_error(peptideIntensityPlot(MSnSetObj = rawMSnSet,
                                      ProteinID = "P03372", 
                                      ProteinName = "ESR1",
                                      combinedIntensities = 1),
                 regexp = errMsg)
    expect_error(peptideIntensityPlot(MSnSetObj = rawMSnSet,
                                      ProteinID = "P03372", 
                                      ProteinName = "ESR1",
                                      combinedIntensities = rawMSnSet),
                 regexp = errMsg)
})
test_that("argument checks - selectedSequence", {
    expect_error(peptideIntensityPlot(rawMSnSet, 
                                      combinedIntensities=protMSnSet,
                                      ProteinID="P03372", 
                                      ProteinName= "ESR1",
                                      selectedSequence = 1),
              regexp = "selectedSequence should be either a character or NULL")
  errMsg <- str_c("selectedSequence: The sequence\\(s\\) provided can not be ",
                  "found in the feature data for P03372")
  expect_error(peptideIntensityPlot(rawMSnSet, 
                                    combinedIntensities=protMSnSet,
                                    ProteinID="P03372", 
                                    ProteinName= "ESR1",
                                    selectedSequence = "[K].1234.[L]"),
               regexp = errMsg)
})
test_that("argument checks - selectedModifications", {
  expect_error(peptideIntensityPlot(rawMSnSet, 
                                    combinedIntensities=protMSnSet,
                                    ProteinID="P03372", 
                                    ProteinName= "ESR1",
                                 selectedSequence = "[K].NVVPLYDLLLEMLDAHR.[L]",
                                    selectedModifications = 1),
          regexp = "selectedModifications should be either a character or NULL")
  errMsg <- str_c("selectedModifications: the modification\\(s\\) provided ",
                  "can not be found in the feature data for P03372")
  expect_error(peptideIntensityPlot(rawMSnSet, 
                                    combinedIntensities=protMSnSet,
                                    ProteinID="P03372", 
                                    ProteinName= "ESR1",
                                 selectedSequence = "[K].NVVPLYDLLLEMLDAHR.[L]",
                                    selectedModifications = "Wibble [N-Term]"),
               regexp = errMsg)
})