# Performs scaling normalization on the intensities within group 

groupScaling <- function(MSnSetObj, scalingFunction=median, 
                         groupingColumn="SampleGroup") {
    checkArg_groupScaling(MSnSetObj, scalingFunction, groupingColumn)
    
    exprs(MSnSetObj) <- as.data.frame(exprs(MSnSetObj)) %>%
        rownames_to_column("PeptideID") %>%
        gather("SampleName", "RawIntensity", -PeptideID) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        rename_at(vars(groupingColumn), ~ "Grouping_column") %>%
        group_by(SampleName) %>%
        mutate(scaledIntensity = scalingFunction(RawIntensity) %>% log()) %>%
        group_by(Grouping_column) %>%
        mutate(meanscaledIntensity = mean(scaledIntensity)) %>%
        ungroup() %>%
        mutate(scalingFactors = exp(scaledIntensity - meanscaledIntensity)) %>%
        mutate(normalizedIntensities = RawIntensity / scalingFactors) %>%
        select(PeptideID, SampleName, normalizedIntensities) %>%
        spread(SampleName, normalizedIntensities) %>%
        arrange(factor(PeptideID, levels = rownames(MSnSetObj))) %>%
        as.data.frame() %>%
        column_to_rownames("PeptideID") %>%
        select(colnames(MSnSetObj)) %>%
        as.matrix()
    return(MSnSetObj)
}

