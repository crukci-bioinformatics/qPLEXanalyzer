# Argument check function
checkArg_coveragePlot <- function(MSnSetObj,
                                  ProteinID,
                                  ProteinName,
                                  fastaFile,
                                  myCol){
    assert_that(is_MSnSet(MSnSetObj), is_PeptideSet(MSnSetObj))
    assert_that("Sequences" %in% colnames(fData(MSnSetObj)),
                msg= 'MSnSetObj feature data must include a "Sequences" column')
    assert_that(is_validProteinId(ProteinID, MSnSetObj, allowNull=FALSE))
    assert_that(is.string(ProteinName))
    assert_that(is.readable(fastaFile))    
    assert_that(length(myCol)==1)
    assert_that(is_validColour(myCol))
}

#' Plot peptide sequence coverage
#' 
#' Computes and displays peptide sequence coverage in proteomics experiment
#' 
#' In the qPLEX-RIME experiment it is imperative for bait protein to have good
#' sequence coverage.  This function plots the protein sequence coverage of the
#' bait protein in the qPLEX-RIME experiment. It requires the fasta sequence
#' file of bait protein as input to generate the plot.
#' 
#' @param MSnSetObj MSnSet: an object of class MSnSet
#' @param ProteinID character: Uniprot ID of the protein
#' @param ProteinName character: name of the protein
#' @param fastaFile character: fasta file of protein sequence
#' @param myCol colour: colour for plotting
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
#' mySequenceFile <- system.file('extdata', 
#'                               "P03372.fasta", 
#'                               package="qPLEXanalyzer")
#' coveragePlot(MSnSet_data, 
#'              ProteinID="P03372", 
#'              ProteinName="ERa", 
#'              fastaFile=mySequenceFile)
#' 
#' @import ggplot2
#' @importFrom Biobase fData
#' @importFrom Biostrings readAAStringSet vmatchPattern
#' @importFrom dplyr distinct filter select pull
#' @importFrom IRanges IRanges reduce start end width 
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @importFrom stringr str_c str_replace_all
#'
#' @export coveragePlot
coveragePlot <- function(MSnSetObj, ProteinID, ProteinName, fastaFile, 
                         myCol="brown") {
    checkArg_coveragePlot(MSnSetObj, ProteinID, ProteinName, fastaFile, 
                             myCol)

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
        pull(Sequences) %>%
        str_replace_all("^\\[.\\]\\.([A-Z]+)\\.\\[.\\]$", "\\1") %>%
        map_dfr(getPosition) %>% 
        distinct()
    
    # get percent coverage
    protWidth <- IRanges::width(Protein_seq)
    coverage <- IRanges(features$start, features$end) %>%
        IRanges::reduce() %>%
        IRanges::width() %>%
        sum()
    Perct <- round(coverage / protWidth * 100, 2)
    SubTitle <- str_c("Number of Unique Peptides: ", 
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

