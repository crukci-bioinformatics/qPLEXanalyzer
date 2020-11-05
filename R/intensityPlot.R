# Argument check function
checkArg_intensityPlot <- function(MSnSetObj,
                                   sampleColours,
                                   title,
                                   colourBy,
                                   transform,
                                   xlab,
                                   trFunc){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is_validSampleColours(sampleColours, colourBy, MSnSetObj))
    assert_that(is.string(title))
    assert_that(is.string(colourBy))
    assert_that(is_validMetadataColumn(colourBy, MSnSetObj))
    assert_that(is.flag(transform))
    assert_that(is.string(xlab))
    assert_that(is.function(trFunc))
    assert_that(length(trFunc(10))==1 & is.numeric(trFunc(10)),
                msg = str_c("trFunc: the specified function should tranform a ",
                            "numeric value into another single numeric value,",
                            "e.g. log2 or sqrt"))
}

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
#' @return An object created by \code{ggplot}
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
#' intensityPlot(MSnSet_data, 
#'               title = "qPLEX_RIME_ER", 
#'               sampleColours = customCols)
#' 
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @importFrom dplyr across everything filter left_join mutate
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect one_of
#'
#' @export intensityPlot
intensityPlot <- function(MSnSetObj, sampleColours=NULL, title="", 
                          colourBy="SampleGroup", transform=TRUE, 
                          xlab="log2(intensity)", trFunc=log2xplus1) {
    if (!transform) { trFunc <- c }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    checkArg_intensityPlot(MSnSetObj, sampleColours, title, colourBy,
                           transform, xlab, trFunc)

    colourBy <- sym(colourBy)

    plotDat <- exprs(MSnSetObj) %>%
        as.data.frame() %>%
        pivot_longer(names_to = "SampleName", 
                     values_to = "Intensity", 
                     everything()) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(!!colourBy, as.factor)) %>%
        mutate(Intensity = trFunc(Intensity)) %>%
        filter(!is.na(Intensity))
    
    
    ggplot(plotDat, aes(x = Intensity, 
                        group = SampleName, 
                        colour = !!colourBy)) +
        stat_density(geom = "line", position = "identity") +
        scale_colour_manual(values = sampleColours) +
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

