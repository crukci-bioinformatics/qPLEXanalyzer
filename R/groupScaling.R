# Argument check function
checkArg_groupScaling <- function(MSnSetObj, scalingFunction, groupingColumn){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is.string(groupingColumn))
    assert_that(is_validMetadataColumn(groupingColumn, MSnSetObj))
}

#' Normalization by scaling within group
#' 
#' Performs scaling normalization on the intensities within group (median or
#' mean)
#' 
#' In this normalization method the central tendencies (mean or median) of the
#' samples within groups are aligned. The argument "groupingColumn" is used to
#' define separate groups to normalize. The function takes one of the column of
#' pData(data) as the variable for classifying group. The default variable is
#' "SampleGroup".  It is imperative in qPLEX-RIME experiment to define IgG as a
#' separate group and normalize it separately from others. You could add a
#' column into the metadata to define this classification.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param scalingFunction function; median or mean
#' @param groupingColumn character; the feature on which groups would be based;
#' default="SampleGroup"
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
#' MSnset_norm <- groupScaling(MSnSet_data, 
#'                             scalingFunction=median, 
#'                             groupingColumn="SampleGroup")
#' 
#' @importFrom Biobase exprs exprs<- pData
#' @importFrom dplyr across arrange group_by left_join mutate select ungroup
#' pull
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr pivot_longer pivot_wider
#'
#' @export groupScaling
groupScaling <- function(MSnSetObj, scalingFunction=median, 
                         groupingColumn="SampleGroup") {
    checkArg_groupScaling(MSnSetObj, scalingFunction, groupingColumn)
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    
    # check none of the samples is entirely missing in intensities
    checNA <- intensities %>%  
        summarize(across(everything(), ~sum(!is.na(.x)))) %>% 
        min() %>% 
        `==`(0)
    if(checNA){
        stop("One or more of the samples is entirely missing in the intensity ",
             "matrix.")
    }
    
    # warning if NAs in intensity matrix
    if(any(is.na(intensities))){
        warning("There are missing values in the intensity matrix. ",
                "These will be omitted in the calculation of scaling factors.")
    }
    
    # get median/mean for each sample
    scaledIntensities <- intensities %>%
        summarize(across(everything(), scalingFunction, na.rm=TRUE)) %>%
        mutate(across(everything(), log)) %>% 
        pivot_longer(names_to = "SampleName", values_to="sInt", everything())
    
    # get group mean for each sample
    meanScaledIntensities <- scaledIntensities %>% 
        left_join(pData(MSnSetObj), by="SampleName") %>% 
        group_by(across(groupingColumn)) %>% 
        mutate(meanscaledIntensity = mean(sInt)) %>% 
        pull(meanscaledIntensity)
    
    scaledIntensities <- pull(scaledIntensities, sInt)
    scalingFactors <- exp(scaledIntensities - meanScaledIntensities)
    normalizedIntensities <- t(t(intensities) / scalingFactors)
    exprs(MSnSetObj) <- normalizedIntensities
    return(MSnSetObj)
}

