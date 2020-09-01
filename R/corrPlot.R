# Correlation plot
corrPlot <- function(MSnSetObj, addValues=TRUE, title="", 
                     low_cor_colour="#FFFFFF", high_cor_colour="#B90505") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("data has to be of class MSnSet..")
    }
    if (!is.logical(addValues)) {
        stop("addValues has to be either TRUE or FALSE..")
    }
    
    col2Cols <- c(low_cor_colour, high_cor_colour)
    cor(exprs(MSnSetObj)) %>%
        as.data.frame() %>%
        rownames_to_column("X") %>%
        pivot_longer(names_to = "Y", values_to = "Cor", -X) %>%
        mutate(addValues = addValues) %>%
        mutate(CorTxt = ifelse(addValues == TRUE, round(Cor, 3), "")) %>%
        ggplot(aes(x = X, y = Y, fill = Cor)) +
        geom_tile(col = "grey") +
        geom_text(aes(label = CorTxt)) +
        scale_fill_gradientn(colors = col2Cols, breaks = seq(0, 1, 0.2)) +
        labs(x = NULL, y = NULL, title = title) +
        guides(fill = guide_colorbar(barheight = 10)) +
        theme(
            aspect.ratio = 1,
            axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5, 
                                       size = 13),
            axis.text.y = element_text(size = 13),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank()
        )
}

