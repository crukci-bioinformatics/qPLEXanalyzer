# PCA plot
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
        mutate_at(vars(colourBy), funs(as.factor))
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

