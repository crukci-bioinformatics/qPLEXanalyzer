# Performs scaling normalization on the intensities within group 

groupScaling <- function(MSnSetObj, scalingFunction=median, 
                         groupingColumn="SampleGroup") {
    checkArg_groupScaling(MSnSetObj, scalingFunction, groupingColumn)
    
    exprs(MSnSetObj) <- as.data.frame(exprs(MSnSetObj)) %>%
        rownames_to_column("PeptideID") %>%
        pivot_longer(names_to = "SampleName", values_to = "RawIntensity", -PeptideID) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        group_by(SampleName) %>%
        mutate(scaledIntensity = scalingFunction(RawIntensity) %>% log()) %>%
        group_by(across(groupingColumn)) %>%
        mutate(meanscaledIntensity = mean(scaledIntensity)) %>%
        ungroup() %>%
        mutate(scalingFactors = exp(scaledIntensity - meanscaledIntensity)) %>%
        mutate(normalizedIntensities = RawIntensity / scalingFactors) %>%
        select(PeptideID, SampleName, normalizedIntensities) %>%
        pivot_wider(names_from = "SampleName", 
                    values_from = "normalizedIntensities") %>%
        arrange(factor(PeptideID, levels = rownames(MSnSetObj))) %>%
        column_to_rownames("PeptideID") %>%
        select(colnames(MSnSetObj)) %>%
        as.matrix()
    return(MSnSetObj)
}

