# MA or volcano plot


#' MA or Volcano Plot
#' 
#' MA or Volcano plot of differential statistics results
#' 
#' Genes determined as significant according to the Log Fold Change and False
#' Discovery Rate cutoffs are highlighted in red.
#' 
#' A user specified selection of genes can be highlighted by passing a character
#' vector of Accessions to the \code{selectedGenes} argument. The contents of
#' this vector will be matched with the Accessions column of the
#' \code{diffstats} object to identify rows to highlight. These will be plotted
#' in blue and labeled with the contents of the \code{GeneSymbol} column. Note
#' that if the \code{GeneSymbol} for a selected gene is missing, no label will
#' be apparent.
#' 
#' @param diffstats list; output of computeDiffStats function
#' @param contrast character; contrast of interest to plot differential
#' statistics results
#' @param title character: title for the plot
#' @param controlGroup character; control group such as IgG
#' @param selectedGenes character: a vector defining genes to plot
#' @param fdrCutOff numeric: False Discovery Rate (adj.P.Val) cut off
#' @param lfcCutOff numeric: Log Fold Change (log2FC) cutoff for
#' @param controlLfcCutOff numeric: only plot genes above controlLogFoldChange
#' cutoff
#' @param plotType character: which type of plot to generate: "MA" or "Volcano"
#' @return An object created by \code{ggplot}
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#' indExpData=c(7:16), Sequences=2, Accessions=6)
#' MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
#' MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
#' contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
#' diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
#' maVolPlot(diffstats, contrast = contrasts, plotType="MA")
#' maVolPlot(diffstats, contrast = contrasts, plotType="Volcano")
#' 
#' @import ggplot2
#' @importFrom Biobase fData
#' @importFrom dplyr arrange bind_rows desc mutate
#' @importFrom magrittr %>%
#'
#' @export maVolPlot
maVolPlot <- function(diffstats, contrast, title="", controlGroup = NULL,
                      selectedGenes=NULL, fdrCutOff=0.05,
                      lfcCutOff=1, controlLfcCutOff=1, plotType="MA") {
    # For plotting we will assign the proteins to one of 7 groups:
    # A - selected (user specified in `selectedGenes`) - highlighted blue
    # B - significant - pass cutoffs - highligted red
    # C - non-significant - fail cutoffs - small and grey
    
    if (!plotType %in% c("MA", "Volcano")) {
        stop("plotType should be 'MA' or 'Volcano'..")
    }
    
    if (!all(selectedGenes%in%fData(diffstats$MSnSetObj)$Accessions)) {
        stop("Some of the genes provided in 'selectedGenes' were not found in ",
             "the data table")
    }

    testSignficant <- function(dat) {
        dat$adj.P.Val <= fdrCutOff &
            !is.na(dat$adj.P.Val) &
            abs(dat$log2FC) >= lfcCutOff &
            abs(dat$controlLogFoldChange) >= controlLfcCutOff
    }
    
    daResTab <- suppressMessages(getContrastResults(
        diffstats = diffstats, contrast = contrast,
        controlGroup = controlGroup
    )) %>%
        bind_rows(data.frame(controlLogFoldChange = vector())) %>%
        mutate(controlLogFoldChange = 
                   ifelse(is.na(controlLogFoldChange), 
                          Inf, 
                          controlLogFoldChange)) %>%
        mutate(group = 
                   ifelse(testSignficant(.), 
                          "Significant", 
                          "Non-significant")) %>%
        mutate(group = ifelse(Accessions %in% selectedGenes, 
                              "Selected", 
                              group)) %>%
        mutate(group = factor(group, 
                              levels = c("Selected", 
                                         "Significant", 
                                         "Non-significant"))) %>%
        arrange(desc(group)) %>%
        mutate(phredPval = -log10(adj.P.Val))
    
    if (plotType == "MA") {
        xFactor <- "AvgIntensity"
        yFactor <- "log2FC"
        xLab <- "average log2(Intensity)"
        yLab <- "log2(Fold Change)"
    }
    if (plotType == "Volcano") {
        xFactor <- "log2FC"
        yFactor <- "phredPval"
        xLab <- "log2(Fold Change)"
        yLab <- "-log10(Adjusted P value)"
    }
    
    xNudge <- diff(range(daResTab[, xFactor])) / 100
    
    ggplot(daResTab, aes_string(
        x = xFactor, 
        y = yFactor, 
        colour = "group",
        size = "group",
        shape = "group",
        alpha = "group",
        fill = "group")) +
        geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
        geom_point() +
        geom_text(
            data = subset(daResTab, group == "Selected"), 
            aes(label = GeneSymbol), 
            hjust = 0,
            size = 3.5,
            nudge_x = xNudge
        ) +
        scale_colour_manual(values = c("black", "black", "gray50"), 
                            drop = FALSE) +
        scale_size_manual(values = c(1.8, 1.5, 0.9), drop = FALSE) +
        scale_shape_manual(values = c(21, 21, 20), drop = FALSE) +
        scale_alpha_manual(values = c(1, 1, 0.6), drop = FALSE) +
        scale_fill_manual(
            values = c("cyan", "red"), limits = levels(daResTab$group)[1:2],
            breaks = levels(droplevels(daResTab$group)), name = ""
        ) +
        labs(x = xLab, y = yLab, title = title) +
        theme_bw() +
        theme(
            text = element_text(size = 16),
            plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        guides(
            colour = "none", size = "none", shape = "none", alpha = "none",
            fill = guide_legend(override.aes = list(shape = 21))
        )
}

