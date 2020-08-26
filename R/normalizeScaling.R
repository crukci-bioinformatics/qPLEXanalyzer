# Performs scaling normalization on the intensities within columns
normalizeScaling <- function(MSnSetObj, scalingFunction, ProteinId = NULL) {
    checkArg_normalizeScaling(MSnSetObj, scalingFunction, ProteinId)
    
    intensities <- as.data.frame(exprs(MSnSetObj))
    intensitiesForScaling <- intensities
    
    if (!is.null(ProteinId)) {
        featuredata <- fData(MSnSetObj)
        ### use protein identifier here
        ind <- which(featuredata$Accessions %in% ProteinId)
        if (!length(ind)) {
            stop("Protein not found")
        }
        intensitiesForScaling <- intensities[ind, ]
    }
    
    scaledIntensities <- intensitiesForScaling %>%
        summarize_all(funs(scalingFunction)) %>%
        mutate_all(funs(log)) %>%
        as.numeric()
    
    scalingFactors <- exp(scaledIntensities - mean(scaledIntensities))
    normalizedIntensities <- t(t(intensities) / scalingFactors)
    exprs(MSnSetObj) <- normalizedIntensities
    return(MSnSetObj)
}

