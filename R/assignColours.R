# Argument check function
# Argument check function
checkArg_assignColours <- function(MSnSetObj, colourBy){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is.string(colourBy))
    assert_that(is_validMetadataColumn(colourBy, MSnSetObj))
}

#' Assigns colours to samples in groups
#' 
#' Assigns colours to samples in groups for plotting
#' 
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param colourBy character: column name from pData(MSnSetObj) to use for
#' coloring samples
#' @return A character vector of colors for samples.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
#' metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#' indExpData=c(7:16), Sequences=2, Accessions=6)
#' sampleColours <- assignColours(MSnSet_data)
#' 
#' @import RColorBrewer
#' @importFrom Biobase pData
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export assignColours
assignColours <- function(MSnSetObj, colourBy="SampleGroup") {
    checkArg_assignColours(MSnSetObj, colourBy)
    
    colourGroups <- as.character(pData(MSnSetObj)[, colourBy]) %>%
        sort() %>% 
        unique()
    len <- length(colourGroups)
    clen <- c(len, 3) %>% max() %>% c(8) %>% min()
    coloursF <- brewer.pal(clen, "Dark2") %>% colorRampPalette()
    sampleColours <- setNames(coloursF(len), colourGroups)
    return(sampleColours)
}

