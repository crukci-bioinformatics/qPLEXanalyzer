#' Hierarchical clustering plot
#' 
#' Computes and displays hierarchical clustering plot for samples in MSnSet
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet
#' @param sampleColours character: a named vector of colors for samples, names
#' should be values of \code{colourBy} column
#' @param colourBy character: column name from \code{pData(MSnSetObj)} to use 
#' for coloring samples
#' @param horizontal logical: define orientation of the dendrogram
#' @param title character: the main title for the dendrogram
#' @return An object created by \code{ggplot}
#' @examples
#' 
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
#' hierarchicalPlot(MSnSet_data, title="qPLEX_RIME_ER")
#' 
#' @import ggplot2
#' @importFrom Biobase exprs pData
#' @importFrom dplyr across left_join mutate
#' @importFrom ggdendro dendro_data label
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @importFrom stats as.dendrogram dist hclust
#' @importFrom utils head
#'
#' @export hierarchicalPlot
hierarchicalPlot <- function(MSnSetObj, sampleColours=NULL, 
                             colourBy="SampleGroup", horizontal=TRUE, 
                             title="") {
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    checkArg_hierarchicalPlot(MSnSetObj, sampleColours, colourBy, horizontal,
                              title)

    dendro.dat <- t(log2xplus1(exprs(MSnSetObj))) %>%
        dist(method = "euclidean") %>%
        hclust() %>%
        as.dendrogram() %>%
        dendro_data()

    colourBy <- sym(colourBy)

    labelDat <- dendro.dat$labels %>%
        mutate(SampleName = as.character(label)) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate(across(!!colourBy, as.factor))

    axisBreaks <- pretty(dendro.dat$segments$yend)[-1] %>% head(-1)
    
    if (horizontal) {
        hj <- 0
        ny <- 1
        ang <- 0
    }
    if (!horizontal) {
        hj <- 1
        ny <- -1
        ang <- 90
    }

    hcPlot <- ggplot(dendro.dat$segment) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(
            data = labelDat,
            aes(x = x, y = y, label = SampleName, colour = !!colourBy),
            hjust = hj, nudge_y = ny, angle = ang
        ) +
        guides(colour = FALSE) +
        scale_colour_manual(values = sampleColours) +
        labs(x = NULL, y = "Distance", title = title)
    if (horizontal) {
        hcPlot <- hcPlot +
            scale_y_reverse(expand = c(0.2, 0), breaks = axisBreaks) +
            coord_flip() +
            theme(
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_blank()
            )
    } else {
        hcPlot <- hcPlot +
            scale_y_continuous(expand = c(0.2, 0), breaks = axisBreaks) +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                plot.title = element_text(hjust = 0.5),
                panel.background = element_blank()
            )
    }
    return(hcPlot)
}

