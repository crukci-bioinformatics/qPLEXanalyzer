testSignficant <- function(dat, cutoffs) {
    isSig <- dat$adj.P.Val <= cutoffs$FDR &
        !is.na(dat$adj.P.Val) &
        abs(dat$log2FC) >= cutoffs$LFC

    if(any(str_detect(colnames(dat), "controlLogFoldChange"))){
        isSig <- isSig & abs(dat$controlLogFoldChange) >= cutoffs$cLFC
    }
    return(isSig)
}

# MA or Volcano Plot
# For plotting we will assign the proteins to one of 3 groups:
# A - selected (user specified in `selectedGenes`) - highlighted blue
# B - significant - pass cutoffs - highligted red
# C - non-significant - fail cutoffs - small and grey

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
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
#'                                metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#'                                indExpData=c(7:16), 
#'                                Sequences=2, 
#'                                Accessions=6)
#' MSnset_norm <- groupScaling(MSnSet_data, scalingFunction=median)
#' MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
#' contrasts <- c(tam.24h_vs_vehicle = "tam.24h - vehicle")
#' diffstats <- computeDiffStats(MSnset_Pnorm, contrasts=contrasts)
#' maVolPlot(diffstats, contrast = contrasts, plotType="MA")
#' maVolPlot(diffstats, contrast = contrasts, plotType="Volcano")
#' 
#' @import ggplot2
#' @importFrom Biobase fData
#' @importFrom dplyr arrange bind_rows case_when desc mutate pull
#' @importFrom magrittr %>%
#' @importFrom tidyr replace_na
#'
#' @export maVolPlot
maVolPlot <- function(diffstats, contrast, title="", controlGroup = NULL,
                      selectedGenes=NULL, fdrCutOff=0.05,
                      lfcCutOff=1, controlLfcCutOff=1, plotType="MA") {
    
    if (!plotType %in% c("MA", "Volcano")) {
        stop("plotType should be 'MA' or 'Volcano'..")
    }
    
    if (!all(selectedGenes%in%fData(diffstats$MSnSetObj)$Accessions)) {
        stop("Some of the genes provided in 'selectedGenes' were not found in ",
             "the data table")
    }
    
    # get contrast results
    results <- suppressMessages(
        getContrastResults(diffstats = diffstats,
                           contrast = contrast,
                           controlGroup = controlGroup
                           ))

    # prepare plot data
    cutOffs <- list(FDR = fdrCutOff, LFC = lfcCutOff, cLFC = controlLfcCutOff)
    grpLevels <- c("Selected", "Significant", "Non-significant")
    daResTab <- results  %>%
        mutate(group = case_when(
            Accessions %in% selectedGenes ~ "Selected",
            testSignficant(., cutOffs) ~ "Significant",
            TRUE ~ "Non-significant")) %>%
        mutate(group = factor(group, levels = grpLevels)) %>%
        arrange(desc(group)) %>%
        mutate(phredPval = -log10(adj.P.Val)) %>%
        mutate(SymbolLab = ifelse(group == "Selected", GeneSymbol, ""))
    
    # set axes for plot type
    if (plotType == "MA") {
        xFactor <- sym("AvgIntensity")
        yFactor <- sym("log2FC")
        xLab <- "average log2(Intensity)"
        yLab <- "log2(Fold Change)"
    }
    if (plotType == "Volcano") {
        xFactor <- sym("log2FC")
        yFactor <- sym("phredPval")
        xLab <- "log2(Fold Change)"
        yLab <- "-log10(Adjusted P value)"
    }
    
    # an x nudge factor for the selected gene labels
    xN <- diff(range(pull(daResTab, !!xFactor))) / 100

    # definitions for scales
    sFill <- c(Selected = "cyan", Significant = "red", `Non-significant` = "gray50")
    sColr <- c(Selected = "black", Significant = "black", `Non-significant` = "gray50")
    sSize <- c(Selected = 1.8, Significant = 1.5, `Non-significant` = 0.9)
    sAlph <- c(Selected = 1, Significant = 1, `Non-significant` = 0.6)
    
    # plot
    ggplot(daResTab, aes(x = !!xFactor, y = !!yFactor)) +
        geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
        geom_point(aes(colour = group,
                       size = group,
                       alpha = group,
                       fill = group),
                   shape = 21) +
        geom_text(aes(label = SymbolLab), hjust = 0, size = 3.5, nudge_x = xN) +
        scale_colour_manual(values = sColr) +
        scale_size_manual(values = sSize, name = "") +
        scale_alpha_manual(values = sAlph) +
        scale_fill_manual(values = sFill, name = "") +
        labs(x = xLab, y = yLab, title = title) +
        theme_bw() +
        theme(
            text = element_text(size = 16),
            plot.title = element_text(size = 14, hjust = 0.5)
        ) +
        guides(
            colour = "none", size = "legend", alpha = "none",
            fill = guide_legend(override.aes = list(shape = 21))
        )
}

