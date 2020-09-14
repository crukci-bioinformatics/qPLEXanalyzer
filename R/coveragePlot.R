# Protein coverage plot


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
#' @param myCol character: colour for plotting
#' @return A protein coverage plot of selected protein for quantitative
#' proteomics data.
#' @examples
#' 
#' data(human_anno)
#' data(exp3_OHT_ESR1)
#' MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
#' indExpData=c(7:16), Sequences=2, Accessions=6)
#' mySequenceFile <- system.file('extdata', "P03372.fasta", package="qPLEXanalyzer")
#' coveragePlot(MSnSet_data, ProteinID="P03372", ProteinName="ERa", fastaFile=mySequenceFile)
#' 
#' @export coveragePlot
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

