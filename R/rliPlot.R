# Relative log intensity plot
rliPlot <- function(MSnSetObj, title="", sampleColours=NULL, 
                    colourBy="SampleGroup", omitIgG=TRUE) {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!colourBy %in% colnames(pData(MSnSetObj))) {
        stop("colourBy must a column names in the pData of the MSnSetObj")
    }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    
    # Remove IgG samples is requested
    if (omitIgG) {
        MSnSetObj <- MSnSetObj[, toupper(MSnSetObj$SampleGroup) != "IGG"]
    }
    exprs(MSnSetObj) %>%
        as.data.frame() %>%
        rownames_to_column("RowID") %>%
        pivot_longer(names_to = "SampleName", 
                     values_to = "Intensity", 
                     -RowID) %>%
        mutate(logInt = log2(Intensity)) %>%
        filter(is.finite(logInt)) %>%
        group_by(RowID) %>%
        mutate(medianLogInt = median(logInt)) %>%
        ungroup() %>%
        mutate(RLI = logInt - medianLogInt) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(one_of(colourBy), as.factor)) %>%
        ggplot(aes_string(x = "SampleName", y = "RLI", fill = colourBy)) +
        geom_boxplot(alpha = 0.6) +
        scale_fill_manual(values = sampleColours, 
                          breaks = names(sampleColours)) +
        labs(x = "Sample", y = "RLI", title = title) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        guides(fill = FALSE)
}
