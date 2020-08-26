# Intensity boxplots
intensityBoxplot <- function(MSnSetObj, title="", sampleColours=NULL, 
                             colourBy="SampleGroup") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!colourBy %in% colnames(pData(MSnSetObj))) {
        stop("colourBy must a column names in the pData of the MSnSetObj")
    }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    
    exprs(MSnSetObj) %>%
        as.data.frame() %>%
        gather("SampleName", "Intensity", everything()) %>%
        mutate(logInt = log2(Intensity)) %>%
        filter(is.finite(logInt)) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate_at(vars(colourBy), funs(as.factor)) %>%
        ggplot() +
        geom_boxplot(aes_string(x = "SampleName", 
                                y = "logInt", 
                                fill = colourBy),
                     alpha = 0.6) +
        scale_fill_manual(values = sampleColours, 
                          breaks = names(sampleColours)) +
        labs(x = "Sample", y = "log2(Intensity)", title = title) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        guides(fill = FALSE)
}

