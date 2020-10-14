context("Convert to MSnSet")
library(qPLEXanalyzer)
data(exp2_Xlink)

# There aren't any missing data in these sets, we need to test the removal
exp2_Xlink$intensities[21, 7] <- NA
exp2_Xlink$intensities[27, 9] <- NA
exp2_Xlink$intensities[34, 11] <- NA
exp2_Xlink$intensities[50, 15] <- NA

# We also need to make sure the reordering of the metadata rows is working
# Currently they are the correct order
exp2_Xlink$metadata <- arrange(exp2_Xlink$metadata, SampleGroup)

exp2MSnSet <- convertToMSnset(exp2_Xlink$intensities,
                              metadata = exp2_Xlink$metadata,
                              indExpData = c(7:16), 
                              Sequences = 2, 
                              Accessions = 6)

exp2MSnSetNA <- convertToMSnset(exp2_Xlink$intensities,
                              metadata = exp2_Xlink$metadata,
                              indExpData = c(7:16), 
                              Sequences = 2, 
                              Accessions = 6,
                              rmMissing = FALSE)

exp2IntCol8 <- exp2_Xlink$intensities[,8]
names(exp2IntCol8) <- paste0("peptide_", seq_along(exp2IntCol8))

# test the function

test_that("Intensity matrix correctly imported", {
  expect_equal(nrow(exp2MSnSet), nrow(exp2_Xlink$intensities) - 4)
  expect_equal(nrow(exp2MSnSetNA), nrow(exp2_Xlink$intensities))
  expect_equal(exprs(exp2MSnSetNA)[,2], exp2IntCol8)
  expect_equal(colnames(exprs(exp2MSnSet)), colnames(exp2_Xlink$intensities)[7:16])
})

test_that("Feature data correctly imported", {
  expect_equal(colnames(fData(exp2MSnSet))[2], "Sequences")
  expect_equal(colnames(fData(exp2MSnSet))[6], "Accessions")
})

# test_that("Metadata correctly imported", {
#  This test is wrong as the script reorders the metadata according to the 
# order of the samples in the intensity matrix
#   expect_equal(colnames(exp2MSnSet), exp2_Xlink$metadata$SampleName)
#  This test is not necessary as `pData(obj) <- x` will fail if the metadata 
# does not match the intensity matrix colnames.
#   expect_equal(ncol(pData(exp2MSnSet)), ncol(exp2_Xlink$metadata))
# })

# test the argument checks

test_that("argument checks - metadata", {
  errMsg <- str_c("Metadata must have columns SampleName, SampleGroup, ",
                  "BioRep, and TechRep")
  expect_error(convertToMSnset(ExpObj = exp2_Xlink$intensities,
                               metadata = data.frame(Wibble=1),
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6),
               regexp = errMsg)
})
test_that("argument checks - Sample Data", {
  errMsg <- str_c("The sample names in the ExpObj columns indicated by indExp ",
                  "do not match the sample names in the metadata table")
  testData <- exp2_Xlink$metadata
  testData$SampleName[1] <- "Wibble"
  expect_error(convertToMSnset(ExpObj = exp2_Xlink$intensities,
                               metadata = testData,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6),
               regexp = errMsg)
  errMsg <- "There are non-numeric values in the intensity data provided"
  testData <- exp2_Xlink$intensities
  testData$FA.rep01 <- as.character(testData$FA.rep01)
  expect_error(convertToMSnset(ExpObj = testData,
                               metadata = exp2_Xlink$metadata,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6),
               regexp = errMsg)
})
test_that("argument checks - Sequences column", {
  errMsg <- str_c("Sequences should be count \\(a single positive integer\\) ",
                  "when type is peptide")
  expect_error(convertToMSnset(ExpObj = exp2_Xlink$intensities,
                               metadata =  exp2_Xlink$metadata,
                               indExpData = c(7:16),
                               Sequences = "A",
                               Accessions = 6),
               regexp = errMsg)
  expect_error(convertToMSnset(ExpObj = exp2_Xlink$intensities,
                               metadata =  exp2_Xlink$metadata,
                               indExpData = c(7:16),
                               Sequences = 1:2,
                               Accessions = 6),
               regexp = errMsg)
  expect_error(convertToMSnset(ExpObj = exp2_Xlink$intensities,
                               metadata =  exp2_Xlink$metadata,
                               indExpData = c(7:16),
                               Sequences = 2,
                               Accessions = 6,
                               type = "protein"),
               regexp = "Sequences should to be NULL when type is protein")
})
