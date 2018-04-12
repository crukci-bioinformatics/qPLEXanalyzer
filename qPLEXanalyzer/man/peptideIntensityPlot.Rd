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
  \item{MSnSetObj}{MSnSet; An object of class MSnSet}
  \item{ProteinID}{Character: Uniprot ID of the protein}
  \item{ProteinName}{Character: Name of the protein}
  \item{combinedIntensities}{MSnSet; An object of class MSnSet}
  \item{selectedSequence}{Character: sequence present in the "Sequences" column in fData(MSnSetObj)}
  \item{selectedModifications}{Character: modification present in the "Modifications" column in fData(MSnSetObj)}
}
\value{
An intensity plot for selected protein and associated peptides.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
print(peptideIntensityPlot(MSnSet_data,combinedIntensities=MSnset_P, ProteinID="P03372", ProteinName= "ESR1"))
}
