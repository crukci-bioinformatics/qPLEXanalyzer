\name{assignColours}
\alias{assignColours}
\title{Assigns colours to samples in groups}
\description{
  Assigns colours to samples in groups for plotting
}
\usage{assignColours(MSnSetObj, colourBy="SampleGroup")}
\arguments{
  \item{MSnSetObj}{MSnSet; an object of class MSnSet}
  \item{colourBy}{character: column name from pData(MSnSetObj) to use for 
  coloring samples}
}
\value{A character vector of colors for samples.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,
metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
sampleColours <- assignColours(MSnSet_data)
}
