context("Coverage plot")
library(qPLEXanalyzer)

data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
                               metadata = exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6)
mySequenceFile <- system.file('extdata', "P03372.fasta", package="qPLEXanalyzer")

plt <- coveragePlot(MSnSet_data, 
                    ProteinID="P03372", 
                    ProteinName="ERa", 
                    fastaFile=mySequenceFile)

test_that("coveragePlot works", {
  vdiffr::expect_doppelganger("coverage plot ERa", plt)
})
