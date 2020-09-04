context("Convert to MSnSet")
library(qPLEXanalyzer)
data(exp2_Xlink)

exp2_Xlink$intensities[21, 7] <- NA
exp2_Xlink$intensities[27, 9] <- NA
exp2_Xlink$intensities[34, 11] <- NA
exp2_Xlink$intensities[50, 15] <- NA

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


test_that("Metadata correctly imported", {
  expect_equal(colnames(exp2MSnSet), exp2_Xlink$metadata$SampleName)
  expect_equal(ncol(pData(exp2MSnSet)), ncol(exp2_Xlink$metadata))
})

