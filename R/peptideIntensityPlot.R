# Argument check function
checkArg_peptideIntensityPlot <- function(MSnSetObj,
                                          ProteinID,
                                          ProteinName,
                                          combinedIntensities,
                                          selectedSequence,
                                          selectedModifications) {
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that(is.string(ProteinID))
    assert_that(is_validProteinId(ProteinID, MSnSetObj, allowNull = FALSE))
    assert_that(
        is.null(combinedIntensities) ||
            (
                is_MSnSet(combinedIntensities) &&
                    is_ProteinSet(combinedIntensities)
            ),
        msg = str_c(
            "combinedIntensities should be either NULL or an ",
            "object of class MSnSet containing summarized ",
            "protein data"
        )
    )
    assert_that(are_validSequences(selectedSequence, MSnSetObj, ProteinID))
    assert_that(are_validModifications(selectedModifications,
                                       MSnSetObj,
                                       ProteinID))
}

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
    checkArg_peptideIntensityPlot(MSnSetObj, ProteinID, ProteinName, 
                                  combinedIntensities, selectedSequence,
                                  selectedModifications)
    
    ## get peptide intensities for plotting
    featureDat <- fData(MSnSetObj) %>% rownames_to_column("PeptideID")
    intensities <- exprs(MSnSetObj) %>%
        as.data.frame() %>%
        rownames_to_column("PeptideID") %>%
        pivot_longer(names_to = "SampleName", 
                     values_to = "Intensity", 
                     -PeptideID) %>%
        left_join(featureDat, by = "PeptideID") %>%
        mutate(logIntensity = log2xplus1(Intensity)) %>%
        filter(Accessions == ProteinID)

    if (!nrow(intensities)) {
        warning("No peptides were found for ", ProteinID)
        return(NULL)
    }
    
    intPlot <- ggplot(intensities, aes(x = SampleName, y = logIntensity)) +
        geom_line(aes(group = PeptideID, colour = Sequences), 
                  size = 0.6, 
                  alpha = 0.5,
                  linetype = 2) +
        geom_point(aes(fill = Sequences), shape = 21, size = 2) +
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
    
    ## plot the selected sequences if present
    if(!is.null(selectedSequence)){    
        seqIntensities <- filter(intensities, Sequences %in% selectedSequence)
        if (!is.null(selectedModifications)) {
            seqIntensities %<>% filter(Modifications %in% selectedModifications)
        }
        
        intPlot <- intPlot +
            geom_line(data = seqIntensities, aes(group = PeptideID), 
                      colour = "#6666FF", 
                      size = 1.2,
                      linetype = 6) +
            geom_point(data = seqIntensities,
                       fill = "#0000FF",
                       shape = 21,
                       size = 2.5)
    }
    
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
        
        intPlot <- intPlot +
            geom_line(data = summIntensities, 
                      aes(group = Accessions),
                      colour = "#AAAAAA",
                      size = 1.2,
                      linetype = 6) +
            geom_point(data = summIntensities,
                       fill = "#888888",
                       shape = 21,
                       size = 2.5)
    }
    
    intPlot
}