context("Merge Peptides")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:2000, ]
rawMSnSet <- convertToMSnset(exp3Int,
                             metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                             indExpData = c(7:16),
                             Sequences = 2,
                             Accessions = 6)

MSnSet_phos <- readRDS("convertToMSnset_phospho_msnset.rds")
protMSnSet <- mergePeptides(rawMSnSet, 
                            summarizationFunction = sum, 
                            annotation = human_anno)

phosProtMSnSet <- mergePeptides(MSnSet_phos, 
                                summarizationFunction = sum, 
                                annotation = human_anno, 
                                keepCols = 7:8)
phosProtMSnSetColname <- mergePeptides(MSnSet_phos, 
                                summarizationFunction = sum, 
                                annotation = human_anno, 
                                keepCols = c("Positions in Master Proteins", 
                                             "Modifications in Master Proteins"))

test_that("Merge peptides works", {
  expect_equal_to_reference(protMSnSet, file="mergePeptides_msnset.rds")  
  expect_equal_to_reference(phosProtMSnSet, file="mergePeptides_phos_msnset.rds")
  expect_equal_to_reference(phosProtMSnSetColname, file="mergePeptides_phos_msnset.rds")
})



