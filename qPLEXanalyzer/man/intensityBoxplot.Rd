\name{intensityBoxplot}
\alias{intensityBoxplot}
\title{Intensity Distribution boxplot}
\description{
  Intensity distribution boxplot of all the samples}
\usage{intensityBoxplot(MSnSetObj, title="", sampleColours=NULL, colourBy="SampleGroup")}
\arguments{
   \item{MSnSetObj}{MSnSet; An object of class MSnSet}
   \item{title}{Character; title of the plot}  
   \item{sampleColours}{Character: A vector of colors for samples}
   \item{colourBy}{Character: column name to use for coloring samples from pData(MSnSetObj)}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
print(intensityBoxplot(MSnSet_data, title = "qPLEX_RIME_ER"))
}
