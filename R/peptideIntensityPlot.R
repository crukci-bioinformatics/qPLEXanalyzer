# Peptide intensity plot
## intensities is a data frame containing peptide intensities with columns for
## each sample summarizedIntensities is a data frame containing summarized
## protein-level intensities protein is the protein for which intensities will
## be plotted samples is a list of samples to use in the plot

peptideIntensityPlot <- function(MSnSetObj, ProteinID, ProteinName,
                                 combinedIntensities=NULL, 
                                 selectedSequence=NULL,
                                 selectedModifications=NULL) {
    if (!is(MSnSetObj, "MSnSet")) {
        stop("MSnSetObj has to be of class MSnSet..")
    }
    if (!is.null(combinedIntensities) && !is(combinedIntensities, "MSnSet")) {
        stop("combinedIntensities has to either NULL or object of class MSnSet")
    }
    
    ## get peptide intensities for plotting
    intensities <- exprs(MSnSetObj) %>%
        as.data.frame() %>%
        rownames_to_column("PeptideID") %>%
        gather(SampleName, Intensity, -PeptideID) %>%
        left_join(fData(MSnSetObj) %>% 
                      rownames_to_column("PeptideID"), "PeptideID") %>%
        mutate(logIntensity = log2xplus1(Intensity)) %>%
        filter(Accessions == ProteinID)
    
    ## get intensities for selected sequences if present
    seqIntensities <- filter(intensities, Sequences %in% selectedSequence)
    if (!is.null(selectedModifications)) {
        seqIntensities %<>% filter(Modifications %in% selectedModifications)
    }
    
    ## get intensity for protein level (summarised) data if present
    summIntensities <- data.frame(SampleName = vector(), 
                                  logIntensity = vector(), 
                                  Accessions = vector())
    if (!is.null(combinedIntensities)) {
        summIntensities <- exprs(combinedIntensities) %>%
            as.data.frame() %>%
            rownames_to_column("Accessions") %>%
            gather(SampleName, Intensity, -Accessions) %>%
            left_join(fData(combinedIntensities), "Accessions") %>%
            mutate(logIntensity = log2xplus1(Intensity)) %>%
            filter(Accessions == ProteinID)
    }
    
    if (!nrow(intensities)) {
        warning("No peptides were found for ", ProteinID)
        return(NULL)
    }
    
    ggplot(intensities, aes(x = SampleName, y = logIntensity)) +
        geom_line(aes(group = PeptideID, colour = Sequences), 
                  size = 0.6, 
                  alpha = 0.5,
                  linetype = 2) +
        geom_point(aes(fill = Sequences), shape = 21, size = 2) +
        ## plot the selected sequences if present
        geom_line(data = seqIntensities, aes(group = PeptideID), 
                  colour = "#6666FF", 
                  size = 1.2,
                  linetype = 6) +
        geom_point(data = seqIntensities,
                   fill = "#0000FF",
                   shape = 21,
                   size = 2.5) +
        ## plot the summarised intensities if present
        geom_line(data = summIntensities, aes(group = Accessions),
                  colour = "#AAAAAA",
                  size = 1.2,
                  linetype = 6) +
        geom_point(data = summIntensities,
                   fill = "#888888",
                   shape = 21,
                   size = 2.5) +
        labs(y = "log2(Intensity)", title = ProteinName) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            legend.position = "none"
        )
}

