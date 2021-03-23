#### Argument checking functions ###############################################
## Custom assertions ###########################################################

#' @import assertthat
#' @importFrom Biobase fData pData
#' @importFrom magrittr %>%
#' @importFrom purrr map_lgl
#' @importFrom stringr str_c str_detect

## Checks on raw data ##########################################################

# check metadata
is_validMetadata <- function(metadata){
    columns <- c("SampleName","SampleGroup","BioRep","TechRep")
    all(columns%in%colnames(metadata))
}
on_failure(is_validMetadata) <- function(call, env) {
    "Metadata must have columns SampleName, SampleGroup, BioRep, and TechRep"
}

# check the sample intensity data ####
is_validSampleData <- function(ExpObj, metadata, indExpData){
    sampleData <- ExpObj[, indExpData]
    assert_that(all(colnames(sampleData)%in%metadata$SampleName) &
                    all(metadata$SampleName%in%colnames(sampleData)),
                msg = str_c("The sample names in the ExpObj columns indicated",
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
                    msg = str_c("Sequences should be count (a single positive ",
                                "integer) when type is peptide"))
    }else{
        assert_that(is.null(Sequences),
                    msg = "Sequences should to be NULL when type is protein")
    }
}

# check the annotation table
is_validAnnotationData <- function(annotation){
    assert_that(is.data.frame(annotation))
    columns <- c("Accessions", "Gene", "Description", "GeneSymbol")
    all(columns%in%colnames(annotation))
}
on_failure(is_validAnnotationData) <- function(call, env) {
    paste0("annotation must have the columns ",
           "Accessions, Gene, Description and GeneSymbol")
}

## Checks on the MSnSet object #################################################

# check MSnSetObj
is_MSnSet <- function(MSnSetObj){
    is(MSnSetObj, "MSnSet")
}
on_failure(is_MSnSet) <- function(call, env){
    "MSnSetObj has to be of class MSnSet"
}

# check MSnSet Object is a peptide level data set ####
is_PeptideSet <- function(MSnSetObj){
    has_name(fData(MSnSetObj), "Sequences")
}
on_failure(is_PeptideSet) <- function(call, env){
    "This MSnSet is not a peptide data set"
}

# check the provided protein ID is in the MSnSetObj
is_validProteinId <- function(ProteinID, MSnSetObj, allowNull=TRUE){
    if(allowNull){
        assert_that(is.character(ProteinID) | is.null(ProteinID),
                    msg=paste0("ProteinID should be a character or NULL"))
        is.null(ProteinID) || ProteinID%in%fData(MSnSetObj)$Accessions
    }else{
        assert_that(is.string(ProteinID))
        ProteinID%in%fData(MSnSetObj)$Accessions
    }
}
on_failure(is_validProteinId) <- function(call, env){
    "The ProteinID provided is not found the MSnset feature data"
}

# check a metadata column given in `varName` exists
is_validMetadataColumn <- function(metaColumn, MSnSetObj){
    metaColumn%in%colnames(pData(MSnSetObj))
}
on_failure(is_validMetadataColumn) <- function(call, env){
    varNam <- deparse(call$metaColumn)
    colNam <- eval(call$metaColumn, env)
    str_c(varNam, ": column ", colNam, " not found in the MSnset metadata")
}

# check the column indices are valid
is_validColNumber <- function(cols, tab){
    max(cols) <= ncol(tab)
}

# check the column names are valid
is_validColName <- function(cols, tab){
    all(cols%in%colnames(tab))
}

# check if the columns are in the fData
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

# check selected peptide sequences are in the fData for the specified protein
are_validSequences <- function(Sequence, MSnSetObj, ProteinID){
    assert_that(is.character(Sequence) | is.null(Sequence),
                msg="selectedSequence should be either a character or NULL")
    subMSnSet <- MSnSetObj[fData(MSnSetObj)$Accessions==ProteinID,]
    is.null(Sequence) || Sequence%in%fData(subMSnSet)$Sequences
}
on_failure(are_validSequences) <- function(call, env){
    varNam <- deparse(call$Sequence)
    protAcc <- eval(call$ProteinID, env)
    str_c(varNam, 
          ": The sequence(s) provided can not be found in the feature data for ",
          protAcc)
}

# check a selected modification is in the fData for the specified protein
are_validModifications <- function(Modifications, MSnSetObj, ProteinID){
    assert_that(is.character(Modifications) | is.null(Modifications),
                msg="selectedModifications should be either a character or NULL")
    subMSnSet <- MSnSetObj[fData(MSnSetObj)$Accessions==ProteinID,]
    is.null(Modifications) || all(Modifications%in%fData(subMSnSet)$Modifications)
}
on_failure(are_validModifications) <- function(call, env){
    varNam <- deparse(call$Modification)
    protAcc <- eval(call$ProteinID, env)
    str_c(varNam, 
          ": the modification(s) provided can not be found in the feature data ",
          "for ", protAcc)
}

# check the control column index
is_validControlColumn <- function(controlInd, MSnSetObj){
    assert_that(is.numeric(controlInd)|is.null(controlInd),
                msg = "controlInd has to be either numeric or NULL")
    is.null(controlInd) || all(controlInd <= ncol(MSnSetObj))
}
on_failure(is_validControlColumn) <- function(call, env){
    "controlInd includes indexes for columns not present in the MSnSet"
}

# check the MSnSet is a protein level data set
is_ProteinSet <- function(MSnSetObj){
    !str_detect(rownames(MSnSetObj)[1], "peptide")
}
on_failure(is_ProteinSet) <- function(call, env){
    "This MSnSet is not a summarized protein data set"
}

# check the batch effect column
is_validBatchEffect <- function(batchEffect, MSnSetObj){
    assert_that(is.character(batchEffect) | is.null(batchEffect),
                msg = "batchEffect has to be either character or NULL")
    is.null(batchEffect) || all(batchEffect%in%colnames(pData(MSnSetObj)))
}
on_failure(is_validBatchEffect) <- function(call, env){
    "batchEffect includes columns not found the MSnset metadata"
}

## Checks for scaling functions ################################################

# Check scaling function is appropriate
is_validScalingFunction <- function(scalingFunction){
    assert_that(is.function(scalingFunction))
    identical(scalingFunction, mean) |
        identical(scalingFunction, BiocGenerics::mean) |
        identical(scalingFunction, median)
}
on_failure(is_validScalingFunction) <- function(call, env){
    "scalingFunction should be mean or median'"
}

## Checks for IRS function ##########################################

# check the IRSname exists in Samplegroup ####
is_validIRSname <- function(IRSname, MSnSetObj){
  assert_that(is.string(IRSname))
  IRSname%in%pData(MSnSetObj)$SampleGroup
}
on_failure(is_validIRSname) <- function(call, env){
  "The IRSname provided is not found the MSnset Samplegroup"
}

## Checks for summarisation functions ##########################################

# check the summarisation function is appropriate
is_validSummarizationFunction <- function(summarizationFunction){
    assert_that(is.function(summarizationFunction))
    length(summarizationFunction(seq(10))) == 1
}
on_failure(is_validSummarizationFunction) <- function(call, env){
    "summarizationFunction should be a summary function, e.g. mean or sum'"
}

## Checks on diffstats object ##################################################

# check list is a computeDiffStats output
is_validDiffstats <- function(diffstats){
    c("MSnSetObj", "fittedLM", "fittedContrasts") %>% 
        map_lgl(has_name, x=diffstats) %>% 
        all()
}
on_failure(is_validDiffstats) <- function(call, env){
    "diffstats is not a valid output of the computeDiffStats function"
}

# check contrast requested is valid
is_validContrast <- function(contrast, diffstats){
    assert_that(is.string(contrast))
    contrasts_available <- colnames(diffstats$fittedContrasts$coefficients)
    contrast%in%contrasts_available
}
on_failure(is_validContrast) <- function(call, env){
    contr <- eval(call$contrast, env)
    diffstats <- eval(call$diffstats, env)
    contr_avail <- diffstats$fittedContrasts$coefficients %>%
      colnames() %>%
      str_c(collapse="\n       ")
    str_c("'", contr, "' is not a valid contrast. Available contrasts:\n",
          "       ", contr_avail)
}

# check the control group
is_validControlGroup <- function(controlGroup, diffstats){
    assert_that(is.string(controlGroup) | is.null(controlGroup),
                msg = paste0("controlGroup has to be a string ", 
                             "(character vector of length 1) or NULL"))
    is.null(controlGroup) ||
        controlGroup%in%colnames(diffstats$fittedLM$coefficients)
}
on_failure(is_validControlGroup) <- function(call, env){
  argNam <- call$controlGroup
  ctrlGrp <- eval(argNam, env)
    str_c(argNam, ": '", ctrlGrp, "'",
           " is not found in the diffstats object. It should be one ", 
           "of the SampleGroups in the original MSnSet object.")
}

## Checks for plotting functions ###############################################

# check string is a valid colour

is_validColour <- function(colour){
    tests <- map_lgl(colour, ~tryCatch(is.matrix(col2rgb(.x)),
                                       error = function(e) FALSE))
    all(tests)
}
on_failure(is_validColour) <- function(call, env) {
    varNam <- deparse(call$colour)
    colour <- eval(call$colour, env)
    tests <- map_lgl(colour, ~tryCatch(is.matrix(col2rgb(.x)),
                                       error = function(e) FALSE))
    str_c(
        str_c(varNam, ": ", colour[!tests], " is not a valid colour."),
        collapse="\n       ")
}

# check sample colours for plotting
is_validSampleColours <- function(sampleColours, colourBy, MSnSetObj){
    assert_that(is_validColour(sampleColours))
    colnm <- names(sampleColours)
    colbynm <- pData(MSnSetObj)[[colourBy]] %>% unique()
    assert_that(!is.null(colnm), 
                msg = str_c("sampleColours should be a named vector, where the ",
                            "names are values in the ", colourBy, " column of ",
                            "the MSnSet object"))
    assert_that(length(colbynm)==length(colnm),
                msg = str_c("sampleColours should provide a colour for each ",
                            "value in the ", colourBy, " column of the MSnSet ",
                            "object"))
    assert_that(all(colbynm%in%colnm) && all(colnm%in%colbynm),
                msg = str_c("The names of sampleColours do not match the values ",
                            "in the ", colourBy, " column of the MSnSet object"))
}
