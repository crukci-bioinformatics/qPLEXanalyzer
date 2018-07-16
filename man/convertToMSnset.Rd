\name{convertToMSnset}
\alias{convertToMSnset}
\title{Converts proteomics TMT intensity data to MSnSet}
\description{
  Converts processed TMT peptide intensities to MSnSet
}
\usage{convertToMSnset(ExpObj, metadata, indExpData, Sequences, Accessions, type="peptide", rmMissing=TRUE)}
\arguments{
  \item{ExpObj}{data.frame; a data.frame consisting of quantitative peptide intensities and peptide annotation}
  \item{metadata}{data.frame; a data.frame describing the samples}
  \item{indExpData}{numeric; a numeric vector indicating the column indexes of intensities in ExpObj}
  \item{Sequences}{numeric; a numeric value indicating the index of column consisting of peptide sequence in ExpObj}
  \item{Accessions}{numeric; a numeric value indicating the index of column consisting of protein accession in ExpObj}
  \item{type}{character; what type of data set to create, either 'peptide' or 'protein'}
  \item{rmMissing}{logical; TRUE or FALSE to indicate whether to remove missing data or not}
}
\details{
This function builds an object of class MSnSet from a dataframe consisting of 
quantitative proteomics intensities data and a meta-data describing the samples 
information. This function creates an MSnSet object from the intensities and metadata file. 
The metadata must contain "SampleName", "SampleGroup", "BioRep" and "TechRep" columns.
The function can be used for either peptide intensities or data that has already been summarized to protein level. The type argument should be set to 'protein' for the latter.
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
}
