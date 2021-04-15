context("Internal Reference Scaling normalisation")
library(qPLEXanalyzer)

data(human_anno)
data(ER_ARID1A_KO_MCF7)

processMSnBatch <- function(batchName, annot=human_anno, dat=ER_ARID1A_KO_MCF7){
    intNam <- str_c("intensities_", batchName)
    metaNam <- str_c("metadata_", batchName)
    intMax <- ncol(dat[[intNam]])-1
    convertToMSnset(dat[[intNam]], 
                    dat[[metaNam]],
                    indExpData=c(7:intMax),
                    Sequences=2,
                    Accessions=6) %>% 
        normalizeScaling(scalingFunction = median) %>% 
        summarizeIntensities(summarizationFunction = sum, 
                             annotation = annot) %>% 
        updateFvarLabels(label=batchName) %>% 
        updateSampleNames(label=batchName)
}

MSnset_SET1 <- processMSnBatch("Set1")
MSnset_SET2 <- processMSnBatch("Set2")
MSnset_SET3 <- processMSnBatch("Set3")
MSnset_SET4 <- processMSnBatch("Set4")
MSnset_SET5 <- processMSnBatch("Set5")

MSnset_comb <- MSnbase::combine(MSnset_SET1,
                                MSnset_SET2,
                                MSnset_SET3,
                                MSnset_SET4,
                                MSnset_SET5)
tokeep <- which(complete.cases(fData(MSnset_comb)))
MSnset_comb <- MSnset_comb[tokeep,]
sampleNames(MSnset_comb) <- pData(MSnset_comb)$SampleName
fData(MSnset_comb) <- fData(MSnset_comb)[,c(2,3,6)]
colnames(fData(MSnset_comb)) <- c("Sequences","Modifications", "Accessions")
MSnset_comb_corr <- IRSnorm(MSnset_comb, IRSname="Ref", groupingColumn="Run")

# Test the function

# For the main object we'll check each of the exprs, fData and pData slots

test_that("Internal reference scaling works - exprs", {
    expect_equal_to_reference(exprs(MSnset_comb_corr), 
                              file="IRSnorm_exprs.rds")
})

test_that("Internal reference scaling works - fData", {
    expect_equal_to_reference(fData(MSnset_comb_corr), 
                              file="IRSnorm_fdata.rds")
})

test_that("Internal reference scaling works - pData", {
    expect_equal_to_reference(pData(MSnset_comb_corr), 
                              file="IRSnorm_pdata.rds")
})

# Check that it works if the batch column is a character vector
# We only need to check the exprs slot

MSnset_comb$Run <- str_c("Plex", MSnset_comb$Run)
MSnset_comb_corr <- IRSnorm(MSnset_comb, IRSname="Ref", groupingColumn="Run")

test_that("Internal reference scaling works - exprs", {
    expect_equal_to_reference(exprs(MSnset_comb_corr), 
                              file="IRSnorm_exprs.rds")
})

# Test the argument checks

test_that("argument checks - MSnset", {
    expect_error(IRSnorm(1), regexp = "MSnSetObj has to be of class MSnSet")
})

test_that("argument checks - IRSname", {
    expect_error(IRSnorm(MSnset_comb, IRSname="Wibble"),
            regexp = "The IRSname provided is not found the MSnset Samplegroup")
})

test_that("argument checks - groupingColumn", {
    expect_error(IRSnorm(MSnset_comb, IRSname="Ref", groupingColumn=1),
                 regexp = "groupingColumn is not a string")
    expect_error(IRSnorm(MSnset_comb, IRSname="Ref", groupingColumn="Wibble"),
      regexp = "groupingColumn: column Wibble not found in the MSnset metadata")
})
