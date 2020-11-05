context("Coverage plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
rawMSnSet <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)
mySequenceFile <- system.file('extdata', 
                              "P03372.fasta", 
                              package="qPLEXanalyzer")

plt <- coveragePlot(rawMSnSet, 
                    ProteinID="P03372", 
                    ProteinName="ERa", 
                    fastaFile=mySequenceFile)

# test the function

test_that("coveragePlot works", {
  vdiffr::expect_doppelganger("coverage plot ERa", plt)
})

# test the argument checks

test_that("argument checks - MSnset", {
    expect_error(coveragePlot(MSnSetObj = 1, 
                              ProteinID="P03372", 
                              ProteinName="ERa", 
                              fastaFile=mySequenceFile),
                 regexp = "MSnSetObj has to be of class MSnSet")
    data(human_anno)
    protMSnSet <- summarizeIntensities(rawMSnSet, sum, human_anno)
    expect_error(coveragePlot(MSnSetObj = protMSnSet, 
                              ProteinID="P03372", 
                              ProteinName="ERa", 
                              fastaFile=mySequenceFile),
                 regexp = "This MSnSet is not a peptide data set")
})
test_that("argument checks - ProteinId", {
    expect_error(coveragePlot(MSnSetObj = rawMSnSet, 
                              ProteinID = NULL, 
                              ProteinName="ERa", 
                              fastaFile=mySequenceFile),
                 regexp = "ProteinID is not a string")
    expect_error(coveragePlot(MSnSetObj = rawMSnSet, 
                              ProteinID = "Wibble", 
                              ProteinName="ERa", 
                              fastaFile=mySequenceFile),
         regexp = "The ProteinID provided is not found the MSnset feature data")
})
test_that("argument checks - myCol", {
    expect_error(coveragePlot(MSnSetObj = rawMSnSet, 
                                           ProteinID="P03372", 
                                           ProteinName="ERa", 
                                           fastaFile=mySequenceFile,
                                           myCol = "wibble"), 
                 regexp = "myCol: wibble is not a valid colour.")
})
