# Protein coverage plot
coveragePlot <- function(MSnSetObj, ProteinID, ProteinName, fastaFile, 
                         myCol="brown") {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!is(ProteinID, "character")) {
        stop("ProteinID has to be of class character..")
    }
    if (!is(ProteinName, "character")) {
        stop("ProteinName has to be of class character..")
    }
    if (!is(fastaFile, "character")) {
        stop("fastaFile has to be of class character..")
    }
    if (!is(myCol, "character")) {
        stop("myCol has to be of class character..")
    }
    if (!"Sequences" %in% colnames(fData(MSnSetObj))) {
        stop('MSnSetObj feature data must include a "Sequences" column')
    }
    
    ## read protein sequence from fastafile
    Protein_seq <- readAAStringSet(fastaFile)
    
    ## extract peptide sequence from MsnSet object and match peptide sequence 
    ## with protein sequence and store co-ordinates
    getPosition <- function(peptideSeq, ProteinSeq=Protein_seq) {
        vmatchPattern(peptideSeq, ProteinSeq) %>%
            as.data.frame() %>%
            select(start, end) %>%
            return()
    }
    features <- fData(MSnSetObj) %>%
        filter(Accessions == ProteinID) %>%
        use_series(Sequence) %>%
        gsub("^\\[.\\]\\.([A-Z]+)\\.\\[.\\]$", "\\1", .) %>%
        map_dfr(getPosition)
    
    features <- unique(features)
    
    # get percent coverage
    protWidth <- width(Protein_seq)
    coverage <- GRanges("feature", IRanges(features$start, features$end)) %>%
        reduce() %>%
        width() %>%
        sum()
    Perct <- round(coverage / protWidth * 100, 2)
    SubTitle <- paste0("Number of Unique Peptides: ", 
                       nrow(features), 
                       "\n% Coverage: ", 
                       Perct)
    
    # set the tick positions for the plot
    nTicks <- min(c(7, ceiling(protWidth / 50) + 1))
    brkTicks <- round(seq(0, protWidth, length.out = nTicks), 0)
    
    features %>%
        distinct() %>%
        ggplot() +
        geom_rect(aes(xmin = start - 1, xmax = end, ymin = 0, ymax = 10), 
                  fill = myCol) +
        geom_rect(xmin = 0,
                  xmax = protWidth,
                  ymin = 0,
                  ymax = 10,
                  colour = "black",
                  fill = NA, 
                  size = 0.15) +
        labs(title = ProteinName, subtitle = SubTitle) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),
            plot.subtitle = element_text(hjust = 0.5),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_x_continuous(limits = c(0, protWidth), breaks = brkTicks) +
        scale_y_continuous(limits = c(0, 10), 
                           breaks = c(0, 10), 
                           expand = c(0, 0))
}

