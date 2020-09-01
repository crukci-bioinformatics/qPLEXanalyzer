context("Peptide intensity plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")
MSnset_P <- readRDS("summarizeIntensities_msnset.rds")

# standard
plt1 <- peptideIntensityPlot(MSnSet_data, 
                     combinedIntensities=MSnset_P, 
                     ProteinID="P03372", 
                     ProteinName= "ESR1")

# no combined intensity
plt2 <- peptideIntensityPlot(MSnSet_data, 
                           ProteinID="P03372", 
                           ProteinName= "ESR1")

# selected sequence
plt3 <- peptideIntensityPlot(MSnSet_data, 
                             combinedIntensities=MSnset_P,
                             ProteinID="P03372", 
                             ProteinName= "ESR1",
                             selectedSequence = "[K].GMEHLYSMK.[C]")

# selected sequence & modification
plt4 <- peptideIntensityPlot(MSnSet_data, 
                             combinedIntensities=MSnset_P,
                             ProteinID="P03372", 
                             ProteinName= "ESR1",
                             selectedSequence = "[K].GMEHLYSMK.[C]",
                             selectedModifications = "2xTMT6plex [N-Term; K9]")


test_that("Peptide plot works", {
  vdiffr::expect_doppelganger("Standard", plt1)
  vdiffr::expect_doppelganger("No combined intensity", plt2)
  vdiffr::expect_doppelganger("With selected sequence", plt3)
  vdiffr::expect_doppelganger("With selected modification", plt4)
})
