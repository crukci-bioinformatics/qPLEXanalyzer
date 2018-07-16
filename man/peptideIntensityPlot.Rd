\name{peptideIntensityPlot}
\alias{peptideIntensityPlot}
\title{Plot peptide intensities}
\description{
  Plots all the peptide intensities for the selected protein
}
\usage{peptideIntensityPlot(MSnSetObj, ProteinID, ProteinName, 
                                 combinedIntensities=NULL, selectedSequence=NULL, 
                                 selectedModifications=NULL)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet containing peptide level intensities}
  \item{ProteinID}{character: Uniprot ID of the protein}
  \item{ProteinName}{character: name of the protein}
  \item{combinedIntensities}{MSnSet; an object of class MSnSet containing protein level intensities}
  \item{selectedSequence}{character: sequence present in the "Sequences" column in fData(MSnSetObj)}
  \item{selectedModifications}{character: modification present in the "Modifications" column in fData(MSnSetObj)}
}
\details{
  Providing a summarised protein level MSnSet object to the 'combinedIntensities' argument will add a summed protein
  intensity trace to the plot along with the peptide intensities.
}
\value{
An intensity plot for selected protein and associated peptides.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
peptideIntensityPlot(MSnSet_data, combinedIntensities=MSnset_P, ProteinID="P03372", ProteinName= "ESR1")
}
