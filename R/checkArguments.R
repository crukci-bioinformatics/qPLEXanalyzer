#### Argument checking functions ###############################################

## Custom assertions ###########################################################

# check metadata ####
is_validMetadata <- function(metadata){
  assert_that(is.data.frame(metadata))
  columns <- c("SampleName","SampleGroup","BioRep","TechRep")
  all(columns%in%colnames(metadata))
}
on_failure(is_validMetadata) <- function(call, env) {
  "Metadata must include the columns SampleName, SampleGroup, BioRep, and TechRep"
}

# check MSnSetObj ####
is_MSnSet <- function(MSnSetObj){
  class(MSnSetObj) == "MSnSet"
}
on_failure(is_MSnSet) <- function(call, env){
  "MSnSetObj has to be of class MSnSet"
}

# check the annotation table ####
is_validAnnotationData <- function(annotation){
    assert_that(is.data.frame(annotation))
    columns <- c("Protein", "Gene", "Description", "GeneSymbol")
    all(columns%in%colnames(annotation))
}
on_failure(is_validAnnotationData) <- function(call, env) {
    "annotation must have the columns Protein, Gene, Description and GeneSymbol"
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
    are_equal(scalingFunction, mean) | are_equal(scalingFunction, median)
}
on_failure(is_validScalingFunction) <- function(call, env){
    "scalingFunction should be mean or median'"
}

# check the provided protein ID is in the MSnSetObj ####
is_validProteinId <- function(ProteinID, MSnSetObj){
    assert_that(is.string(ProteinID))
    ProteinID%in%fData(MSnSetObj)$Accessions
}
on_failure(is_validProteinId) <- function(call, env){
    "The ProteinID provided is not found the MSnset feature data"
}

# check the grouping column exists ####
is_validGroupingColumn <- function(groupingColumn, MSnSetObj){
  assert_that(is.string(groupingColumn))
  groupingColumn%in%colnames(pData(MSnSetObj))
}
on_failure(is_validGroupingColumn) <- function(call, env){
  "The grouping column provided is not found the MSnset metadata"
}

# check the summarisation function is appropriate ####
is_validSummarizationFunction <- function(summarizationFunction){
  assert_that(is.function(summarizationFunction))
  length(summarizationFunction(1:10)) == 1
}
on_failure(is_validSummarizationFunction) <- function(call, env){
  "summarizationFunction should be a summary function, e.g. mean or sum'"
}

# check the control column index ####
is_validControlColumn <- function(controlInd, MSnSetObj){
  assert_that(is.numeric(controlInd)|is.null(controlInd),
              msg = "Sequences has to be either numeric or NULL")
  is.null(controlInd) || all(controlInd <= ncol(MSnSetObj))
}
on_failure(is_validControlColumn) <- function(call, env){
  "controlInd includes indexes for columns not present in the MSnSet"
}

# check the MSnSeis a protein level data set ####
is_ProteinSet <- function(MSnSetObj){
    !grepl("peptide", rownames(MSnSetObj))
}
on_failure(is_PeptideSet) <- function(call, env){
    "This MSnSet is not a summarized protein data set"
}

# check the batch effect column ####
is_validBatchEffect <- function(batchEffect, MSnSetObj){
  assert_that(is.character(batchEffect))
  all(batchEffect%in%colnames(pData(MSnSetObj)))
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
  contrasts <- colnames(diffstats$fittedContrasts$coefficients)
  contrast%in%contrasts
}
on_failure(is_validContrast) <- function(call, env){
    str_c("'", call$contrast, "' is not a valid contrast")
}

# check the control group #####
is_validControlGroup <- function(controlGroup, diffstats){
  assert_that(is.string(controlGroup))
  controlGroup%in%colnames(diffstats$fittedLM$coefficients)
}
on_failure(is_validControlGroup) <- function(call, env){
  str_c(call$controlGroup, 
        " is not found in the diffstats object. ", 
        "It should be one of the SampleGroups in the original MSnSet object.")
}

################################################################################

checkArg_convertToMSnset <- function(ExpObj, metadata, indExpData, Sequences, 
                                     Accessions, type, rmMissing){
    assert_that(is.data.frame(ExpObj))
    assert_that(is_validMetadata(metadata))
    assert_that(is.numeric(indExpData), noNA(indExpData))
    assert_that(is.count(Sequences)|is.null(Sequences),
                msg = "Sequences has to be either numeric or NULL")
    assert_that(is.count(Accessions))
    assert_that(type%in%c("peptide", "protein"), 
                msg="type must be either 'peptide' or 'protein'")
    assert_that(is.flag(rmMissing))
}

checkArg_summarizeIntensities <- function(MSnSetObj, summarizationFunction, 
                                          annotation){
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that(is_validSummarizationFunction(summarizationFunction))
    assert_that(is_validAnnotationData(annotation))
}


checkArg_normalizeQuantiles <- function(MSnSetObj){
    assert_that(is_MSnSet(MSnSetObj))
}

checkArg_normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

checkArg_groupScaling <- function(MSnSetObj, scalingFunction, groupingColumn){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is_validGroupingColumn(groupingColumn, MSnSetObj))
}

checkArg_rowScaling <- function(MSnSetObj, scalingFunction){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
}

checkArg_regressIntensity <- function(MSnSetObj, controlInd, ProteinId){
    assert_that(is_MSnSet(MSnSetObj), is_ProteinSet(MSnSetObj))
    assert_that(is_validControlColumn(controlInd, MSnSetObj))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

checkArg_computeDiffStats <- function(MSnSetObj, batchEffect, transform, 
                                      contrasts, trend, robust){
    assert_that(is_MSnSet(MSnSetObj), is_ProteinSet(MSnSetObj))
    assert_that(is_validBatchEffect(batchEffect, MSnSetObj))
    assert_that(is.flag(transform))
    assert_that(is.character(contrasts))
    assert_that(is.flag(trend))
    assert_that(is.flag(robust))
}

checkArg_getContrastResults <- function(diffstats, contrast, controlGroup, 
                                        transform, writeFile){
    assert_that(is_validDiffStats(diffstats))
    assert_that(is_validContrast(checkContrast))
    assert_that(is_validControlGroup(checkControlGroup))
    assert_that(is.flag(transform))
    assert_that(is.flag(writeFile))
}


################################################################################
