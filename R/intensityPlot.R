# Intensity distribution plot
## intensities is a data frame containing columns for each sample
## sampleColours is a named vector that maps samples to colours
intensityPlot <- function(MSnSetObj, sampleColours=NULL, title="", 
                          colourBy="SampleGroup", transform=TRUE, 
                          xlab="log2(intensity)", trFunc=log2xplus1) {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!transform) {
        trFunc <- c
    }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    exprs(MSnSetObj) %>%
        as.data.frame() %>%
        pivot_longer(names_to = "SampleName", 
                     values_to = "Intensity", 
                     everything()) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(one_of(colourBy), as.factor)) %>%
        mutate(Intensity = trFunc(Intensity)) %>%
        filter(!is.na(Intensity)) %>%
        ggplot(aes_string(x = "Intensity", 
                          group = "SampleName", 
                          colour = colourBy)) +
        stat_density(geom = "line", position = "identity") +
        scale_colour_manual(values = sampleColours, 
                            breaks = names(sampleColours)) +
        labs(x = xlab, title = title) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "cm"),
            legend.text = element_text(size = 11),
            legend.justification = c(1, 1),
            legend.position = c(1, 1),
            legend.background = element_rect(fill = "transparent")
        )
}

