\name{assignColours}
\alias{assignColours}
\title{Assigns colours to samples in groups}
\description{
  Assigns colours to samples in groups for plotting
}
\usage{assignColours(data, groupColourSchemes)}
\arguments{
  \item{data}{MSnSet; An object of class MSnSet}
  \item{groupColourSchemes}{character: character vector of colors for samples}
}
\value{A character vector of colors for samples.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
sampleColours <- assignColours(MSnset_data, groupColourSchemes = c("Blues","Greens","Purples",
"Reds","YlOrBr"))
}
