# Argument check function
checkArg_normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validScalingFunction(scalingFunction))
    assert_that(is_validProteinId(ProteinId, MSnSetObj))
}

#' Normalization by scaling
#' 
#' Performs scaling normalization on the peptide/protein intensities (median or
#' mean)
#' 
#' In this normalization method the central tendencies (mean or median) of the
#' samples are aligned.  The central tendency for each sample is computed and
#' log transformed. A scaling factor is determined by subtracting from each
#' central tendency the mean of all the central tendencies.  The raw intensities
#' are then divided by the scaling factor to get normalized intensities.
#' 
#' The intensities can also be normalized based on the peptide intensities of
#' a selected protein. For this the argument "ProteinId" allows you to define
#' the protein that will be used for scaling the intensities.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param scalingFunction function; median or mean
#' @param ProteinId character; protein Id
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
#' MSnset_norm <- normalizeScaling(MSnSet_data, scalingFunction=median)
#' 
#' @importFrom Biobase exprs exprs<- fData
#' @importFrom dplyr across everything mutate summarize
#' @importFrom magrittr %>%
#'
#' @export normalizeScaling
normalizeScaling <- function(MSnSetObj, 
                             scalingFunction=median, 
                             ProteinId = NULL) {
    checkArg_normalizeScaling(MSnSetObj, scalingFunction, ProteinId)
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    intensitiesForScaling <- intensities

    # do we want to normalize against a specific protein
    protNorm <- !is.null(ProteinId)
    if (protNorm) {
        featuredata <- fData(MSnSetObj)
        ### use protein identifier here
        ind <- which(featuredata$Accessions %in% ProteinId)
        if (!length(ind)) {
            stop("Protein not found")
        }
        intensitiesForScaling <- intensities[ind, ]
    }
    
    # check none of the samples is entirely missing in intensitiesForScaling
    checNA <- intensitiesForScaling %>%  
        summarize(across(everything(), ~sum(!is.na(.x)))) %>% 
        min() %>% 
        `==`(0)
    if(checNA & !protNorm){
        stop("One or more of the samples is entirely missing in the intensity ",
             "matrix.")
    }
    if(checNA & protNorm){
        stop("One or more of the samples is entirely missing in the intensity ",
             "matrix for the protein(s) provided.")
    }
    
    # warning if NAs in intensity matrix
    if(any(is.na(intensitiesForScaling))){
        warning("There are missing values in the intensity matrix. ",
                "These will be omitted in the calculation of scaling factors.")
    }
    
    scaledIntensities <- intensitiesForScaling %>%
        summarize(across(everything(), scalingFunction, na.rm = TRUE)) %>%
        mutate(across(everything(), log)) %>%
        as.numeric()
    
    scalingFactors <- exp(scaledIntensities - mean(scaledIntensities))
    normalizedIntensities <- t(t(intensities) / scalingFactors)
    exprs(MSnSetObj) <- normalizedIntensities
    return(MSnSetObj)
}

