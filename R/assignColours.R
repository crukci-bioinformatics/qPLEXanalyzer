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
    if (len < 3) {
        sampleColours <- setNames(brewer.pal(3, "Dark2")[seq_len(len)], 
                                  colourGroups)
    } else if (len <= 8) {
        sampleColours <- setNames(brewer.pal(len, "Dark2"), colourGroups)
    } else {
        coloursF <- brewer.pal(8, "Dark2") %>% colorRampPalette()
        sampleColours <- setNames(coloursF(len), colourGroups)
    }
    return(sampleColours)
}

