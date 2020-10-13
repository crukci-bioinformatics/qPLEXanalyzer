# Argument check function
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

#' Converts proteomics TMT intensity data to MSnSet
#' 
#' Converts processed TMT peptide intensities to MSnSet
#' 
#' This function builds an object of class MSnSet from a dataframe consisting
#' of quantitative proteomics intensities data and a meta-data describing the
#' samples information. This function creates an MSnSet object from the
#' intensities and metadata file.  The metadata must contain "SampleName",
#' "SampleGroup", "BioRep" and "TechRep" columns. The function can be used for
#' either peptide intensities or data that has already been summarized to
#' protein level. The type argument should be set to 'protein' for the latter.
#' 
#' @param ExpObj data.frame; a data.frame consisting of quantitative peptide
#' intensities and peptide annotation
#' @param metadata data.frame; a data.frame describing the samples
#' @param indExpData numeric; a numeric vector indicating the column indexes of
#' intensities in ExpObj
#' @param Sequences numeric; a numeric value indicating the index of column
#' consisting of peptide sequence in ExpObj
#' @param Accessions numeric; a numeric value indicating the index of column
#' consisting of protein accession in ExpObj
#' @param type character; what type of data set to create, either 'peptide' or
#' 'protein'
#' @param rmMissing logical; TRUE or FALSE to indicate whether to remove
#' missing data or not
#' @return An object of class \code{MSnSet} (see \code{\link{MSnSet-class}}) 
#' object).
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#'
#' @import MSnbase
#' @importFrom Biobase fData featureNames featureNames<- pData pData<- 
#' sampleNames
#' @importFrom dplyr left_join mutate
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom tibble column_to_rownames tibble
#' @importFrom tidyr drop_na
#'
#' @export convertToMSnset
convertToMSnset <- function(ExpObj, metadata, indExpData, Sequences=NULL, 
                            Accessions, type="peptide", rmMissing=TRUE) {
    checkArg_convertToMSnset(ExpObj, metadata, indExpData, Sequences, 
                             Accessions, type, rmMissing)
    
    colnames(ExpObj)[Accessions] <- "Accessions"
    colnames(ExpObj)[Sequences] <- "Sequences"
    
    if (rmMissing) {
        ExpObj <- ExpObj %>% 
            drop_na(indExpData)
    }
    obj <- readMSnSet2(ExpObj, ecol = indExpData)

    pData(obj) <- tibble(SampleName=sampleNames(obj)) %>% 
        left_join(metadata, by="SampleName") %>% 
        mutate(rowname = SampleName) %>% 
        column_to_rownames()
            
    if (type == "protein") {
        featureNames(obj) <- fData(obj)$Accessions
    } else {
        featureNames(obj) <- str_c("peptide_", featureNames(obj))
    }
    return(obj)
}

