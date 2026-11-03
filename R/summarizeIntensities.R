# Argument check function
checkArg_summarizeIntensities <- function(MSnSetObj, 
                                          summarizationFunction, 
                                          annotation){
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that(is_validSummarizationFunction(summarizationFunction))
    assert_that(is_validAnnotationData(annotation))
}

##### Summarization function ##########

# Functions for summarizing intensity data

# Summarizes multiple peptide measurements for a protein.
# Assumes that there are columns for each of the samples specified and uses
# Accessions column for grouping peptide-level measurements.
# Missing values are not filtered: they are passed to the summarization
# function as they are.
# Typical summarization functions are sum, mean and median.
# The annotation must have the columns "Accessions", "Gene", "Description" and
# "GeneSymbol".
# The MSnSetObj must have column "Sequences" denoting its a peptide dataset



#' Summarizes peptides intensities to proteins
#'
#' Summarizes multiple peptides intensities measurements to protein level.
#'
#' The peptide intensities of each protein (grouped by the "Accessions" column
#' of fData(MSnSetObj)) are combined in each sample with
#' \code{summarizationFunction}. Missing values are passed to this function as
#' they are: with \code{sum}, \code{mean} or \code{median} a protein gets a
#' missing value in a sample if any of its peptides is missing in that sample.
#' The resulting feature data contain the annotation columns and a "Count"
#' column giving the number of distinct peptide sequences (values of the
#' "Sequences" column) of each protein. Proteins whose Accession is not found
#' in the annotation are kept, with missing annotation.
#'
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param summarizationFunction function; method used to aggregate the peptides
#' into proteins. Sum, mean or median
#' @param annotation data.frame; a data.frame of protein annotation of four
#' columns: "Accessions", "Gene", "Description" and "GeneSymbol"
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}})
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
#' 
#' @import MSnbase
#' @importFrom Biobase exprs fData featureNames featureNames<- pData pData<- 
#' sampleNames
#' @importFrom dplyr across bind_cols everything group_by left_join mutate
#' @importFrom dplyr n_distinct select summarize
#' @importFrom tidyselect where
#'
#' @export summarizeIntensities
summarizeIntensities <- function(MSnSetObj, summarizationFunction, annotation) {
    checkArg_summarizeIntensities(MSnSetObj, summarizationFunction, annotation)
    
   summarizedIntensities <- fData(MSnSetObj) |>
       select(Accessions, Sequences) |>
       mutate(across(c(Accessions, Sequences), as.character)) |> 
       bind_cols(as.data.frame(exprs(MSnSetObj))) |> 
       group_by(Accessions) |>
       summarize(across(where(is.numeric), summarizationFunction), 
                 Count=n_distinct(Sequences)) |> 
       left_join(annotation, by = "Accessions") |>
       select(Accessions, colnames(annotation), Count, everything())
    
    obj <- readMSnSet2(summarizedIntensities, ecol = sampleNames(MSnSetObj))
    pData(obj) <- pData(MSnSetObj)
    featureNames(obj) <- fData(obj)$Accessions
    sampleNames(obj) <- pData(obj)$SampleName
    return(obj)
}

