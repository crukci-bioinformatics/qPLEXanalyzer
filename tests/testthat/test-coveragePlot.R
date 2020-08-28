context("Coverage plot")
library(qPLEXanalyzer)

MSnSet_data <- readRDS("convertToMSnset_oht_esr1_msnset.rds")
mySequenceFile <- system.file('extdata', "P03372.fasta", package="qPLEXanalyzer")

plt <- coveragePlot(MSnSet_data, 
                    ProteinID="P03372", 
                    ProteinName="ERa", 
                    fastaFile=mySequenceFile)

test_that("coveragePlot works", {
  vdiffr::expect_doppelganger("coverage plot ERa", plt)
})