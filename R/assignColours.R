# Assigns colours
## check this for selecting colors http://colorbrewer2.org/

assignColours <- function(MSnSetObj, colourBy="SampleGroup") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    colourGroups <- as.character(pData(MSnSetObj)[, colourBy]) %>%
        sort() %>% 
        unique()
    len <- length(colourGroups)
    clen <- c(len, 3) %>% max() %>% c(8) %>% min()
    coloursF <- brewer.pal(clen, "Dark2") %>% colorRampPalette()
    sampleColours <- setNames(coloursF(len), colourGroups)
    return(sampleColours)
}

