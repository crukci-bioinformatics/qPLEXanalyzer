# MA or volcano plot
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

