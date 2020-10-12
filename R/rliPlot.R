#' Relative log intensity plot
#' 
#' Relative log intensity (RLI) plots of the samples within MSnset
#' 
#' An RLI-plot is a boxplot that can be used to visualise unwanted variation in
#' a data set. It is similar to the relative log expression plot developed for
#' microarray analysis - see Gandolfo and Speed (2018).  Rather than examining
#' gene expression, the RLI plot uses the MS intensities for each peptide or
#' the summarised protein intensities.
#' 
#' The column provided to the \code{colourBy} argument will be used to colour
#' the samples. The colours will be determined using the function
#' \code{\link{assignColours}}, alternatively the user may specify a named
#' vector of colours using the \code{sampleColours} argument. The names of the
#' \code{sampleColours} vector should match the unique values in the
#' \code{colourBy} column.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param title character: title for the plot
#' @param sampleColours character: a named vector of colours for samples
#' @param colourBy character: column name to use for colouring samples from
#' pData(MSnSetObj)
#' @param omitIgG logical: whether to remove IgG from the RLI plot
#' @return An object created by \code{ggplot}
#' @references Gandolfo LC, Speed TP (2018) RLE plots: Visualizing unwanted
#' variation in high dimensional data. PLoS ONE 13(2): e0191629.
#' \url{https://doi.org/10.1371/journal.pone.0191629}
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' rliPlot(MSnSet_data, title = "qPLEX_RIME_ER")
#' 
#' # custom colours
#' customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
#' names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
#' rliPlot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)
#' 
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @importFrom dplyr across filter group_by left_join mutate ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom stats median
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect one_of
#'
#' @export rliPlot
rliPlot <- function(MSnSetObj, title="", sampleColours=NULL, 
                    colourBy="SampleGroup", omitIgG=TRUE) {
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    checkArg_rliPlot(MSnSetObj, title, sampleColours, colourBy, omitIgG)

    # Remove IgG samples is requested
    if (omitIgG) {
        MSnSetObj <- MSnSetObj[, toupper(MSnSetObj$SampleGroup) != "IGG"]
    }
    
    colourBy <- sym(colourBy)
    plotDat <- exprs(MSnSetObj) %>%
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
        mutate(across(!!colourBy, as.factor))
    
    plotDat %>%
        ggplot(aes(x = SampleName, y = RLI, fill = !!colourBy)) +
        geom_boxplot(alpha = 0.6) +
        scale_fill_manual(values = sampleColours) +
        labs(x = "Sample", y = "RLI", title = title) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        guides(fill = FALSE)
}
