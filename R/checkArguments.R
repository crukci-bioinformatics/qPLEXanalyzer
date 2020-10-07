#### Argument checking functions ###############################################
## Custom assertions ###########################################################

#' @import assertthat
#' @importFrom Biobase fData pData
#' @importFrom magrittr %>%
#' @importFrom purrr map_lgl
#' @importFrom stringr str_c str_detect

# check string is a valid colour ####

is_validColour <- function(colour){
  tryCatch(is.matrix(col2rgb(colour)),
           error = function(e) FALSE)
}
on_failure(is_validColour) <- function(call, env) {
  varNam <- deparse(call$colour)
  colVal <- eval(call$colour, env)
  str_c(varNam, ": ", colVal, " is not a valid colour.")
}

# check metadata ####
is_validMetadata <- function(metadata){
    assert_that(is.data.frame(metadata))
    columns <- c("SampleName","SampleGroup","BioRep","TechRep")
    all(columns%in%colnames(metadata))
}
on_failure(is_validMetadata) <- function(call, env) {
    "Metadata must have columns SampleName, SampleGroup, BioRep, and TechRep"
}

# check the sample dat ####
is_validSampleData <- function(ExpObj, metadata, indExpData){
    sampleData <- ExpObj[, indExpData]
    assert_that(all(colnames(sampleData)%in%metadata$SampleName)&
                    all(metadata$SampleName%in%colnames(sampleData)),
                msg=paste0("The sample names in the ExpObj columns indicated",
                           " by indExp do not match the sample names in the ",
                           "metadata table"))
    all(map_lgl(sampleData, is.numeric))
}
on_failure(is_validSampleData) <- function(call, env) {
    "There are non-numeric values in the intensity data provided"
}

# check the peptide sequences column ####
is_validSequencesColumn <- function(Sequences, type){
    if(type=="peptide"){
        assert_that(is.count(Sequences), 
                    msg=paste0("Sequences should be count (a single positive ", 
                               "integer) when type is peptide"))
    }else{
        is.null(Sequences)
    }
}
on_failure(is_validSequencesColumn) <- function(call, env){
    "Sequences should to be NULL when type is protein"
}

# check MSnSetObj ####
is_MSnSet <- function(MSnSetObj){
    is(MSnSetObj, "MSnSet")
}
on_failure(is_MSnSet) <- function(call, env){
    "MSnSetObj has to be of class MSnSet"
}

# check the annotation table ####
is_validAnnotationData <- function(annotation){
    assert_that(is.data.frame(annotation))
    columns <- c("Accessions", "Gene", "Description", "GeneSymbol")
    all(columns%in%colnames(annotation))
}
on_failure(is_validAnnotationData) <- function(call, env) {
    paste0("annotation must have the columns ",
           "Accessions, Gene, Description and GeneSymbol")
}

# check MSnSetOBj is a peptide level data set ####
is_PeptideSet <- function(MSnSetObj){
    has_name(fData(MSnSetObj), "Sequences")
}
on_failure(is_PeptideSet) <- function(call, env){
    "This MSnSet is not a peptide data set"
}

# Check scaling function is appropriate ####
is_validScalingFunction <- function(scalingFunction){
    assert_that(is.function(scalingFunction))
    are_equal(scalingFunction, BiocGenerics::mean) |
        are_equal(scalingFunction, mean) | 
        are_equal(scalingFunction, median)
}
on_failure(is_validScalingFunction) <- function(call, env){
    "scalingFunction should be mean or median'"
}

# check the provided protein ID is in the MSnSetObj ####
is_validProteinId <- function(ProteinID, MSnSetObj){
    assert_that(is.character(ProteinID) | is.null(ProteinID),
                msg=paste0("ProteinID should be a charater or NULL"))
    is.null(ProteinID) || ProteinID%in%fData(MSnSetObj)$Accessions
}
on_failure(is_validProteinId) <- function(call, env){
    "The ProteinID provided is not found the MSnset feature data"
}

# check a metadata column given in `varName` exists ####
is_validMetadataColumn <- function(metaColumn, MSnSetObj){
    assert_that(is.string(metaColumn))
    metaColumn%in%colnames(pData(MSnSetObj))
}
on_failure(is_validMetadataColumn) <- function(call, env){
    varNam <- deparse(call$metaColumn)
    colNam <- eval(call$metaColumn, env)
    str_c(varNam, ": column ", colNam, " not found in the MSnset metadata")
}

# check the summarisation function is appropriate ####
is_validSummarizationFunction <- function(summarizationFunction){
    assert_that(is.function(summarizationFunction))
    length(summarizationFunction(seq(10))) == 1
}
on_failure(is_validSummarizationFunction) <- function(call, env){
    "summarizationFunction should be a summary function, e.g. mean or sum'"
}

# check the column indices are valid ####
is_validColNumber <- function(cols, tab){
  max(cols) <= ncol(tab)
}

# check the column names are valid ####
is_validColName <- function(cols, tab){
  all(cols%in%colnames(tab))
}

# check if the columns are in the fData ####
is_validfDataColumn <- function(cols, MSnSetObj){
  tab <- fData(MSnSetObj)
  errMsg1 <- str_c("keepCols should be NULL, a character vector of column ",
                   "names or a numeric vector of column indices.")
  errMsg2 <- str_c("One or more of keepCols exceeds the number of columns in ",
                   "fData(MSnSetObj).")
  errMsg3 <- "One or more of keepCols is not found in fData(MSnSetObj)."
  assert_that(is.numeric(cols) | is.character(cols) | is.null(cols),
              msg = errMsg1)
  if(is.numeric(cols)){ 
    assert_that(is_validColNumber(cols, tab), msg = errMsg2) 
    }
  if(is.character(cols)){ 
    assert_that(is_validColName(cols, tab), msg = errMsg3) 
    }
  TRUE
}

# check the control column index ####
is_validControlColumn <- function(controlInd, MSnSetObj){
    assert_that(is.numeric(controlInd)|is.null(controlInd),
                msg = "controlInd has to be either numeric or NULL")
    is.null(controlInd) || all(controlInd <= ncol(MSnSetObj))
}
on_failure(is_validControlColumn) <- function(call, env){
    "controlInd includes indexes for columns not present in the MSnSet"
}

# check the MSnSet is a protein level data set ####
is_ProteinSet <- function(MSnSetObj){
    !str_detect(rownames(MSnSetObj)[1], "peptide")
}
on_failure(is_ProteinSet) <- function(call, env){
    "This MSnSet is not a summarized protein data set"
}

# check the batch effect column ####
is_validBatchEffect <- function(batchEffect, MSnSetObj){
    assert_that(is.character(batchEffect) | is.null(batchEffect),
                msg = "batchEffect has to be either character or NULL")
    is.null(batchEffect) || all(batchEffect%in%colnames(pData(MSnSetObj)))
}
on_failure(is_validBatchEffect) <- function(call, env){
    "batchEffect includes columns not found the MSnset metadata"
}

# check list is a computeDiffStats output ####
is_validDiffstats <- function(diffstats){
    c("MSnSetObj", "fittedLM", "fittedContrasts") %>% 
        map_lgl(has_name, x=diffstats) %>% 
        all()
}
on_failure(is_validDiffstats) <- function(call, env){
    "diffstats is not a valid output of the computeDiffStats function"
}

# check contrast requested is valid ####
is_validContrast <- function(contrast, diffstats){
    assert_that(is.string(contrast))
    contrasts_available <- colnames(diffstats$fittedContrasts$coefficients)
    contrast%in%contrasts_available
}
on_failure(is_validContrast) <- function(call, env){
    str_c("'", call$contrast, "' is not a valid contrast")
}

# check the control group #####
is_validControlGroup <- function(controlGroup, diffstats){
    assert_that(is.string(controlGroup) | is.null(controlGroup),
                msg = paste0("controlGroup has to be a string ", 
                             "(character vector of length 1) or NULL"))
    is.null(controlGroup) ||
        controlGroup%in%colnames(diffstats$fittedLM$coefficients)
}
on_failure(is_validControlGroup) <- function(call, env){
    str_c(call$controlGroup, 
           " is not found in the diffstats object. It should be one ", 
           "of the SampleGroups in the original MSnSet object.")
}

# Main argument check functions ################################################

checkArg_assignColours <- function(MSnSetObj, colourBy){
  assert_that(is_MSnSet(MSnSetObj))
  assert_that(is_validMetadataColumn(colourBy, MSnSetObj))
}

checkArg_computeDiffStats <- function(MSnSetObj, batchEffect, transform, contrasts, trend, robust){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validBatchEffect(batchEffect, MSnSetObj))
    assert_that(is.flag(transform))
    assert_that(is.character(contrasts))
    assert_that(is.flag(trend))
    assert_that(is.flag(robust))
}

checkArg_convertToMSnset <- function(ExpObj, 
                                     metadata, 
                                     indExpData, 
                                     Sequences, 
                                     Accessions, 
                                     type, 
                                     rmMissing){
    assert_that(is.data.frame(ExpObj))
    assert_that(is_validMetadata(metadata))
    assert_that(is.numeric(indExpData), noNA(indExpData))
    assert_that(is_validSampleData(ExpObj, metadata, indExpData))
    assert_that(type%in%c("peptide", "protein"), 
                msg="type must be either 'peptide' or 'protein'")
    assert_that(is_validSequencesColumn(Sequences, type))
    assert_that(is.count(Accessions))
    assert_that(is.flag(rmMissing))
}

checkArg_corrPlot <- function(MSnSetObj,
                              addValues,
                              title,
                              low_cor_colour,
                              high_cor_colour,
                              textsize){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is.flag(addValues))
    assert_that(is.string(title))
    assert_that(length(low_cor_colour)==1)
    assert_that(is_validColour(low_cor_colour))
    assert_that(length(high_cor_colour)==1)
    assert_that(is_validColour(high_cor_colour))
    assert_that(is.number(textsize))
}

checkArg_getContrastResults <- function(diffstats, 
                                        contrast, 
                                        controlGroup, 
                                        transform, 
                                        writeFile){
    assert_that(is_validDiffstats(diffstats))
    assert_that(is_validContrast(contrast, diffstats))
    assert_that(is_validControlGroup(controlGroup, diffstats))
    assert_that(is.flag(transform))
    assert_that(is.flag(writeFile))
}

checkArg_groupScaling <- function(MSnSetObj, scalingFunction, groupingColumn){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is_validMetadataColumn(groupingColumn, MSnSetObj))
}

checkArg_mergePeptides <- function(MSnSetObj, 
                                   summarizationFunction, 
                                   annotation, 
                                   keepCols){
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that(is_validSummarizationFunction(summarizationFunction))
    assert_that(is_validAnnotationData(annotation))
    assert_that(is_validfDataColumn(keepCols, MSnSetObj))
}

checkArg_normalizeQuantiles <- function(MSnSetObj){
    assert_that(is_MSnSet(MSnSetObj))
}

checkArg_normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

checkArg_regressIntensity <- function(MSnSetObj, controlInd, ProteinId){
    assert_that(is_MSnSet(MSnSetObj), is_ProteinSet(MSnSetObj))
    assert_that(is_validControlColumn(controlInd, MSnSetObj))
    assert_that(is.string(ProteinId))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

checkArg_rowScaling <- function(MSnSetObj, scalingFunction){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
}

checkArg_summarizeIntensities <- function(MSnSetObj, 
                                          summarizationFunction, 
                                          annotation){
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that(is_validSummarizationFunction(summarizationFunction))
    assert_that(is_validAnnotationData(annotation))
}

################################################################################
