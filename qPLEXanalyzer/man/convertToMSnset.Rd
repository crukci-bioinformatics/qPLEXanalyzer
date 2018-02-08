\name{convertToMSnset}
\alias{convertToMSnset}
\title{Converts proteomics TMT intensity data to MSnset}
\description{
  Converts processed TMT intensity data to MSnset
}
\usage{convertToMSnset(data,metadata,indExpData,indFData,rmMissing=TRUE)}
\arguments{
  \item{data}{data frame; An object of class data frame}
  \item{metadata}{A dataframe describing the samples (in lines)}
  \item{indExpData}{A numeric vector indicating the column indexes of intensity in data}
  \item{indFData}{A numeric vector indicating the column indexes of features in data}
  \item{rmMissing}{logical; TRUE or FALSE to indicate whether to remove missing data or not}
}
\details{
This function is a wrapper of createMSnset function of builds an object of class MSnSet 
from a single tabulated-like file for quantitative and meta-data and a dataframe for the 
samples description. This function creates a MSnSet object from the intensities and metadata 
file. The metadata must contain "Experiment","Label", "Bio.Rep","Tech.Rep" and "Analyt.Rep" 
columns. "Experiment" and "Label" must be read in as factors.
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
}
