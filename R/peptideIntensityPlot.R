# Peptide intensity plot
## intensities is a data frame containing peptide intensities with columns for
## each sample summarizedIntensities is a data frame containing summarized
## protein-level intensities protein is the protein for which intensities will
## be plotted samples is a list of samples to use in the plot



#' Plot peptide intensities
#' 
#' Plots all the peptide intensities for the selected protein
#' 
#' Providing a summarised protein level MSnSet object to the
#' \code{combinedIntensities} argument will add a summed protein intensity trace
#' to the plot along with the peptide intensities.
#' 
#' @param MSnSetObj MSnSet; an object of class MSnSet containing peptide level
#' intensities
#' @param ProteinID character: Uniprot ID of the protein
#' @param ProteinName character: name of the protein
#' @param combinedIntensities MSnSet; an object of class MSnSet containing
#' protein level intensities
#' @param selectedSequence character: sequence present in the "Sequences"
#' column in \code{fData(MSnSetObj)}
#' @param selectedModifications character: modification present in the
#' "Modifications" column in \code{fData(MSnSetObj)}
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
#' MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
#' peptideIntensityPlot(MSnSet_data, 
#'                      combinedIntensities=MSnset_P, 
#'                      ProteinID="P03372", 
#'                      ProteinName= "ESR1")
#' 
#' @import ggplot2
#' @importFrom Biobase exprs fData
#' @importFrom dplyr filter left_join mutate
#' @importFrom magrittr %<>% %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @export peptideIntensityPlot
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
        pivot_longer(names_to = "SampleName", 
                     values_to = "Intensity", 
                     -PeptideID) %>%
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
            pivot_longer(names_to = "SampleName", 
                         values_to = "Intensity", 
                         -Accessions) %>%
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

