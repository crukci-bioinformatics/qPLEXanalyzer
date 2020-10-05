context("Merge Peptides")
library(qPLEXanalyzer)

data(human_anno)
data(exp3_OHT_ESR1)
exp3Int <- exp3_OHT_ESR1$intensities_qPLEX1[1:1000, ]
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

# test using column names instead of numbers
phosProtMSnSetColname <- mergePeptides(MSnSet_phos, 
                                summarizationFunction = sum, 
                                annotation = human_anno, 
                                keepCols = c("Positions in Master Proteins", 
                                             "Modifications in Master Proteins"))

# # The MSnSet object contains the MSnbase version in the `processingData` slot
# # This will cause the test to fail if the MSnbase version used in the build
# test_that("Merge peptides works", {
#   expect_equal_to_reference(protMSnSet, file="mergePeptides_msnset.rds")  
#   expect_equal_to_reference(phosProtMSnSet, file="mergePeptides_phos_msnset.rds")
#   expect_equal_to_reference(phosProtMSnSetColname, file="mergePeptides_phos_msnset.rds")
# })

# The function creates an entirely new MSnSet obj, we should compare the samples
# features data, and merged intensities of each object

protTestList <- list(Samples = pData(protMSnSet),
                     Features = fData(protMSnSet),
                     MergedIntensitied = exprs(protMSnSet))

# these two should be identical, it was just the method of selecting the columns
# that changed.
phosProtTestList <- list(Samples = pData(phosProtMSnSet),
                         Features = fData(phosProtMSnSet),
                         MergedIntensitied = exprs(phosProtMSnSet))

phosProtColTestList <- list(Samples = pData(phosProtMSnSetColname),
                         Features = fData(phosProtMSnSetColname),
                         MergedIntensitied = exprs(phosProtMSnSetColname))

test_that("Merge peptides works", {
    expect_equal_to_reference(protTestList, 
                              file="mergePeptides.rds")  
    expect_equal_to_reference(phosProtTestList, 
                              file="mergePeptides_phos.rds")
    expect_equal_to_reference(phosProtColTestList, 
                              file="mergePeptides_phos.rds")
})



