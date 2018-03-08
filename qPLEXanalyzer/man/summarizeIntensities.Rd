\name{summarizeIntensities}
\alias{summarizeIntensities}
\title{Summarizes peptides intensities to proteins}
\description{
  Summarizes multiple peptides measurements to proteins.
}
\usage{summarizeIntensities(MSnSetObj, summarizationFunction, annotation)}
\arguments{
  \item{MSnSetObj}{MSnSet; An object of class MSnSet}
  \item{summarizationFunction}{function; Method used to aggregate the peptides 
  into proteins. Sum, mean or median}
  \item{annotation}{data frame; A dataframe of protein annotation of four columns: 
                        "Protein", "Gene", "Description" and "GeneSymbol"}
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
MSnset_P <- summarizeIntensities(MSnSet_data, sum, human_anno)
}
