\name{peptideIntensityPlot}
\alias{peptideIntensityPlot}
\title{Plot peptide intensities}
\description{
  Plots all the peptide intensities for the selected protein
}
\usage{peptideIntensityPlot(peptideIntensities, combinedIntensities=NULL, protein,
                                 title = "", ylab = "log2(intensity)",
                                 minIntensity = NA, maxIntensity = NA,
                                 selectedSequence = NULL, selectedModifications = NULL)}
\arguments{
  \item{peptideIntensities}{MSnSet; An object of class MSnSet}
  \item{combinedIntensities}{MSnSet; An object of class MSnSet}
  \item{protein}{Character: UniprotID of the protein to be plotted}
  \item{title}{Character: title for the plot}
  \item{ylab}{Character: text for y-axis}
  \item{minIntensity}{Add text here}
  \item{maxIntensity}{Add text here}
  \item{selectedSequence}{Add text here}
  \item{selectedModifications}{Add text here}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
MSnset_norm <- groupScaling(MSnset_data,func=median, Grp="Label")
MSnset_Pnorm <- summarizeIntensities(MSnset_norm, sum, human_anno)
gene <- "ESR1"
protein <- "P03372"
print(peptideIntensityPlot(MSnset_norm,combinedIntensities=MSnset_Pnorm, protein, title = gene))
}
