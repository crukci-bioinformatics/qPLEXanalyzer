\name{normalizeScaling}
\alias{normalizeScaling}
\title{Normalization by scaling}
\description{
  Performs scaling normalization on the peptide/protein intensities (median or mean)
}
\usage{normalizeScaling(MSnSetObj,func,Protein = NULL)}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{func}{function; median or mean}
  \item{Protein}{character; protein Id}
}
\details{
  In this normalization method the central tendencies (mean or median) of the samples are aligned. 
  The central tendency for each sample is computed and log transformed. A scaling factor is determined 
  by subtracting from each central tendency the mean of all the central tendencies.  The raw intensities
  are then divided by the scaling factor to get normalized ones. 
  The intensities can also be normalized based on the peptides intensities of a selected protein. For this
  the argument "Protein" allows you to define the protein that will be used for scaling the intensities.
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
MSnset_norm <- normalizeScaling(MSnSet_data,func=median)
}
