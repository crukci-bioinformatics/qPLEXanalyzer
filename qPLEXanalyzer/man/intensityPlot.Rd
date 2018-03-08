\name{intensityPlot}
\alias{intensityPlot}
\title{Intensity Distribution Plot}
\description{
  Intensity distribution plot of all the samples}
\usage{intensityPlot(MSnSetObj, sampleColours=NULL, title="", colourBy="SampleGroup", 
                          transform=TRUE, xlab="log2(intensity)",
                          trFunc=log2xplus1)}
\arguments{
   \item{MSnSetObj}{MSnSet; An object of class MSnSet}
   \item{sampleColours}{Character: A vector of colors for samples}
   \item{title}{Character: Title for the plot}
   \item{colourBy}{Character: column name to use for coloring samples from pData(MSnSetObj)}
   \item{transform}{Logical: Whether to log transform intensities}
   \item{xlab}{Character: Text for x-axis}
   \item{trFunc}{func: internal helper function for log transformation}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
print(intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER"))
}
