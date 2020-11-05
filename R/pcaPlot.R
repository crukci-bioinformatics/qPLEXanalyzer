# Argument check function
checkArg_pcaPlot <- function(MSnSetObj,
                             omitIgG,
                             sampleColours,
                             transFunc,
                             transform,
                             colourBy,
                             title,
                             labelColumn,
                             labelsize,
                             pointsize,
                             x.nudge,
                             x.PC){
    assert_that(is_MSnSet(MSnSetObj))
    assert_that(is.flag(omitIgG))
    assert_that(is_validSampleColours(sampleColours, colourBy, MSnSetObj))
    assert_that(is.function(transFunc))
    assert_that(length(transFunc(10))==1 & is.numeric(transFunc(10)),
                msg = str_c("transFunc: the specified function should tranform ",
                            " a numeric value into another single numeric value,",
                            "e.g. log2 or sqrt"))
    assert_that(is.flag(transform))
    assert_that(is.string(colourBy))
    assert_that(is_validMetadataColumn(colourBy, MSnSetObj))
    assert_that(is.string(title))
    assert_that(is.string(labelColumn))
    assert_that(is_validMetadataColumn(labelColumn, MSnSetObj))
    assert_that(is.number(labelsize))
    assert_that(is.number(pointsize))
    assert_that(is.number(x.nudge))
    assert_that(is.count(x.PC))
}

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
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @importFrom dplyr across left_join mutate
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom stats na.omit prcomp
#' @importFrom stringr str_c
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect one_of
#'
#' @export pcaPlot
pcaPlot <- function(MSnSetObj, omitIgG=FALSE, sampleColours=NULL,
                    transFunc=log2xplus1, transform=TRUE, 
                    colourBy="SampleGroup", title="", labelColumn="BioRep",
                    labelsize=4, pointsize=4, x.nudge=4, x.PC=1) {
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    checkArg_pcaPlot(MSnSetObj, omitIgG, sampleColours, transFunc, transform, 
                     colourBy, title, labelColumn, labelsize, pointsize, 
                     x.nudge, x.PC)
    
    ## Remove IgG samples is requested
    if (omitIgG) {
        MSnSetObj <- MSnSetObj[, toupper(MSnSetObj$SampleGroup) != "IGG"]
    }
    
    intensities <- exprs(MSnSetObj) %>%
        as.data.frame() %>%
        na.omit()
    
    if (!nrow(intensities)) { return(NULL) }
    
    if(transform) {  intensities <- transFunc(intensities) }
    
    pca <- intensities %>% t() %>% prcomp()
    pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
    
    colourBy <- sym(colourBy)
    plotDat <- as.data.frame(pca$x) %>%
        rownames_to_column("SampleName") %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(!!colourBy, as.factor))
    
    xLab <- str_c("PC", x.PC, ", ", pcaVariance[x.PC], "% variance")
    yLab <- str_c("PC", x.PC + 1, ", ", pcaVariance[x.PC + 1], "% variance")
    
    xPC <- str_c("PC", x.PC) %>% sym()
    yPC <- str_c("PC", x.PC + 1) %>% sym()
    
    pcPlot <- ggplot(plotDat, aes(x = !!xPC, y = !!yPC)) +
        geom_point(aes(fill = !!colourBy),
                   pch = 21, 
                   colour = "black",
                   size = pointsize) +
        scale_fill_manual(values = sampleColours) +
        labs(x = xLab,y = yLab, fill = NULL, title = title) +
        theme_bw() +
        theme(
            text = element_text(size = 14),
            plot.title = element_text(size = 14, hjust = 0.5),
            aspect.ratio = 1
        )
    
    if (!is.null(labelColumn)) {
        labelColumn <- sym(labelColumn) 
        pcPlot <- pcPlot +
            geom_text(aes(label = !!labelColumn),
                      hjust = 0,
                      size = labelsize,
                      nudge_x = x.nudge)
    }
    
    pcPlot
}

