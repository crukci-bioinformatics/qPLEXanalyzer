# Intensity distribution plot
## intensities is a data frame containing columns for each sample
## sampleColours is a named vector that maps samples to colours


#' Intensity Distribution Plot
#' 
#' Intensity distribution plot of all the samples
#' 
#' The column provided to the \code{colourBy} argument will be used to colour
#' the samples. The colours will be determined using the function
#' \code{\link{assignColours}}, alternatively the user may specify a named
#' vector of colours using the \code{sampleColours} argument. The names of the
#' \code{sampleColours} vector should match the unique values in the
#' \code{colourBy} column.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param sampleColours character: a vector of colors for samples
#' @param title character: title for the plot
#' @param colourBy character: column name from pData(MSnSetObj) to use for
#' coloring samples
#' @param transform logical: whether to log transform intensities
#' @param xlab character: label for x-axis
#' @param trFunc func: internal helper function for log transformation
#' @return An intensity distribution plot for quantitative proteomics data.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1, 
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER")
#' 
#' # custom colours
#' customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
#' names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
#' intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)
#' 
#' @export intensityPlot
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

