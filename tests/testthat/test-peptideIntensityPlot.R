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

