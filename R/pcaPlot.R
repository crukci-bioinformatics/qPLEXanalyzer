# PCA plot


#' PCA plot
#'
#' PCA plots of the samples within MSnset
#'
#' The column provided to the "colourBy" argument will be used to colour the
#' samples. The colours will be determined using the function
#' \link{assignColours}, alternatively the user may specify a named vector of
#' colours using the "sampleColours" argument. The names of the "sampleColours"
#' vector should match the unique values in the "colourBy" column.
#'
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param omitIgG Logical: whether to remove IgG from the PCA plot
#' @param sampleColours character: A named vector of colours for samples
#' @param transFunc func: internal helper function for log transformation
#' @param transform logical: whether to log transform intensities
#' @param colourBy character: column name to use for colouring samples from
#'   \code{pData(MSnSetObj)}
#' @param title character: title for the plot
#' @param labelColumn character: column name from \code{pData(MSnSetObj)} to use
#'   for labelling points on the plot
#' @param labelsize numeric: size of the labels
#' @param pointsize numeric: size of plotting points
#' @param x.nudge numeric: distance to move labels along the x-axis away from
#'   the plotting points
#' @param x.PC numeric: The principle component to plot on the x-axis; the
#'   following PC will be plotted on the y-axis
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
#' exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
#' pcaPlot(MSnSet_data, omitIgG = TRUE, labelColumn = "BioRep")
#'
#' # custom colours and PC2 v PC3
#' customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
#' names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
#' pcaPlot(MSnSet_data, 
#'         omitIgG = TRUE, 
#'         labelColumn = "BioRep", 
#'         sampleColours = customCols, 
#'         x.PC=2)
#'
#' @export pcaPlot
pcaPlot <- function(MSnSetObj, omitIgG=FALSE, sampleColours=NULL,
                    transFunc=log2xplus1, transform=TRUE, 
                    colourBy="SampleGroup", title="", labelColumn="BioRep",
                    labelsize=4, pointsize=4, x.nudge=4, x.PC=1) {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("data has to be of class MSnSet..")
    }
    if (!colourBy %in% colnames(pData(MSnSetObj))) {
        stop("colourBy must a column names in the pData of the MSnSetObj")
    }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    
    ## Remove IgG samples is requested
    if (omitIgG) {
        MSnSetObj <- MSnSetObj[, toupper(MSnSetObj$SampleGroup) != "IGG"]
    }
    if (!transform) {
        transFunc <- as.data.frame
    }
    intensities <- exprs(MSnSetObj) %>%
        as.data.frame() %>%
        na.omit() %>%
        transFunc()
    if (!nrow(intensities)) {
        return(NULL)
    }
    pca <- intensities %>% t() %>% prcomp()
    pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
    plotDat <- as.data.frame(pca$x) %>%
        rownames_to_column("SampleName") %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(one_of(colourBy), as.factor))
    xPC <- paste0("PC", x.PC)
    yPC <- paste0("PC", x.PC + 1)
    ggplot(plotDat, 
           aes_string(x = xPC, 
                      y = yPC, 
                      fill = colourBy, 
                      label = labelColumn)) +
        geom_point(pch = 21, colour = "black", size = pointsize) + {
            if (!is.null(labelColumn)) {
                geom_text(hjust = 0, size = labelsize, nudge_x = x.nudge)
            }
        } +
        scale_fill_manual(values = sampleColours, 
                          breaks = names(sampleColours)) +
        labs(
            x = paste0(xPC, ", ", pcaVariance[x.PC], "% variance"),
            y = paste0(yPC, ", ", pcaVariance[x.PC + 1], "% variance"),
            fill = NULL, title = title
        ) +
        theme_bw() +
        theme(
            text = element_text(size = 14),
            plot.title = element_text(size = 14, hjust = 0.5),
            aspect.ratio = 1
        )
}

