# Hierachical clustering plot
hierarchicalPlot <- function(MSnSetObj, sampleColours=NULL, 
                             colourBy="SampleGroup", horizontal=TRUE, 
                             title="") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!colourBy %in% colnames(pData(MSnSetObj))) {
        stop("colourBy must a column names in the pData of the MSnSetObj")
    }
    if (is.null(sampleColours)) {
        sampleColours <- assignColours(MSnSetObj, colourBy = colourBy)
    }
    
    dendro.dat <- t(log2xplus1(exprs(MSnSetObj))) %>%
        dist(method = "euclidean") %>%
        hclust() %>%
        as.dendrogram() %>%
        dendro_data()
    labelDat <- dendro.dat$labels %>%
        mutate(SampleName = as.character(label)) %>%
        left_join(pData(MSnSetObj), "SampleName") %>%
        mutate_at(vars(colourBy), funs(as.factor))
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
            aes_string(x = "x", 
                       y = "y", 
                       label = "SampleName", 
                       colour = colourBy),
            hjust = hj, nudge_y = ny, angle = ang
        ) +
        guides(colour = FALSE) +
        scale_fill_manual(values = sampleColours, 
                          breaks = names(sampleColours)) +
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

