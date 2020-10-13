# Argument check function
checkArg_mergePeptides <- function(MSnSetObj, 
                                   summarizationFunction, 
                                   annotation, 
                                   keepCols){
  assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
  assert_that(is_validSummarizationFunction(summarizationFunction))
  assert_that(is_validAnnotationData(annotation))
  assert_that(is_validfDataColumn(keepCols, MSnSetObj))
}

#' Merge identical peptides intensities
#' 
#' Merge identical peptides to single peptide intensity. This function is
#' especially useful for phosphopeptide analysis.
#' 
#' Rows of the intensity matrix with identical peptide sequences are merged by
#' summarising the intensities using \code{summarizationFunction}.
#' 
#' Columns specified with \code{keepCols} are retained in the final output.
#' Non-unique entries in different rows are concatenated with ';'.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param summarizationFunction function; method used to aggregate the
#' peptides. sum, mean or median
#' @param annotation data.frame; a data.frame of protein annotation of four
#' columns: "Accessions", "Gene", "Description" and "GeneSymbol"
#' @param keepCols a vector of additional columns from fData(MSnSetObj) to
#' keep.  either be a numeric vector of column numbers of a character vector of
#' column names
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
#' MSnset_P <- mergePeptides(MSnSet_data, sum, human_anno)
#' 
#' @import MSnbase
#' @importFrom Biobase exprs fData featureNames featureNames<- pData pData<- 
#' sampleNames
#' @importFrom Biobase sampleNames
#' @importFrom dplyr across bind_cols everything group_by left_join mutate n 
#' @importFrom dplyr select summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom tidyselect all_of
#'
#' @export mergePeptides
mergePeptides <- function(MSnSetObj, 
                          summarizationFunction, 
                          annotation, 
                          keepCols=NULL) {
  checkArg_mergePeptides(MSnSetObj, summarizationFunction, annotation, keepCols)
  
  concatUnique <- function(x){ unique(x) %>% str_c(collapse=";") }
  
  summarizedIntensities <- fData(MSnSetObj) %>%
    select(Accessions, Sequences, all_of(keepCols)) %>%
    mutate(across(everything(), as.character)) %>% 
    mutate(phosseqid = str_c(Sequences, "_", Accessions)) %>%  
    bind_cols(as.data.frame(exprs(MSnSetObj))) %>% 
    group_by(phosseqid, Accessions) %>%
    select(-Sequences) %>% 
    summarize(across(where(is.character), concatUnique),
              across(where(is.numeric), summarizationFunction), 
              Count=n())   %>% 
    ungroup() %>% 
    left_join(annotation, by = "Accessions") %>%
    select(Accessions, 
           colnames(annotation), 
           Count, 
           everything())

  obj <- readMSnSet2(summarizedIntensities, ecol = sampleNames(MSnSetObj))
  pData(obj) <- pData(MSnSetObj)
  featureNames(obj) <- fData(obj)$phosseqid
  return(obj)
}

