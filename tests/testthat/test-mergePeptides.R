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

# R CMD check (and devtools::check()) is using a different locale to 
# devtools::test(). This leads to the summarised table being sorted differently
# specifying the locale ensure they are both the same and we don't get errors.
Sys.setlocale(category = "LC_COLLATE", locale = "C")
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

# test the function

test_that("Merge peptides works", {
     expect_equal_to_reference(protTestList,
                              file="mergePeptides.rds") 
    expect_equal_to_reference(phosProtTestList,
                              file="mergePeptides_phos.rds")
    expect_equal_to_reference(phosProtColTestList,
                              file="mergePeptides_phos.rds")
})

# test the argument checks

test_that("argument checks - MSnset", {
    expect_error(mergePeptides(MSnSetObj = 1,
                               summarizationFunction = sum,
                               annotation = human_anno),
                 regexp = "MSnSetObj has to be of class MSnSet")
    protMSnSet <- summarizeIntensities(rawMSnSet, sum, human_anno)
    expect_error(mergePeptides(MSnSetObj = protMSnSet,
                               summarizationFunction = sum,
                               annotation = human_anno),
                 regexp = "This MSnSet is not a peptide data set")
})
test_that("argument checks - summarizationFunction", {
    expect_error(mergePeptides(MSnSetObj = rawMSnSet,
                               summarizationFunction = 1,
                               annotation = human_anno),
                 regexp = "summarizationFunction is not a function")
    errMsg <- str_c("summarizationFunction should be a summary function, ",
                    "e.g. mean or sum'")
    expect_error(mergePeptides(MSnSetObj = rawMSnSet,
                               summarizationFunction = sqrt,
                               annotation = human_anno),
                 regexp = errMsg)
})
test_that("argument checks - annotation", {
    expect_error(mergePeptides(MSnSetObj = rawMSnSet,
                               summarizationFunction = sum,
                               annotation = 1),
                 regexp = "annotation is not a data frame")
    testAnnot <- rename(human_anno, Wibble = Accessions)
    errMsg <- str_c("annotation must have the columns Accessions, Gene, ",
                    "Description and GeneSymbol")
    expect_error(mergePeptides(MSnSetObj = rawMSnSet,
                               summarizationFunction = sum,
                               annotation = testAnnot),
                 regexp = errMsg)
})
test_that("argument checks - keepCols", {
    errMsg <- str_c("keepCols should be NULL, a character vector of column ",
                    "names or a numeric vector of column indices.")
    expect_error(mergePeptides(MSnSetObj = MSnSet_phos,
                               summarizationFunction = sum,
                               annotation = human_anno,
                               keepCols = NA),
                 regexp = errMsg)
    errMsg <- str_c("One or more of keepCols exceeds the number of columns ",
                    "in fData\\(MSnSetObj\\).")
    expect_error(mergePeptides(MSnSetObj = MSnSet_phos,
                               summarizationFunction = sum,
                               annotation = human_anno,
                               keepCols = c(1,2,7,12)),
                 regexp = errMsg)
    errMsg <- "One or more of keepCols is not found in fData\\(MSnSetObj\\)."
    expect_error(mergePeptides(MSnSetObj = MSnSet_phos,
                               summarizationFunction = sum,
                               annotation = human_anno,
                               keepCols = c("Positions in Master Proteins",
                                            "Modifications in Master Proteins",
                                            "Wibble")),
                 regexp = errMsg)
})

