\name{intensityPlot}
\alias{intensityPlot}
\title{Intensity Distribution Plot}
\description{
  Intensity distribution plot of all the samples}
\usage{intensityPlot(MSnSetObj, sampleColours=NULL, title="", colourBy="SampleGroup", 
                          transform=TRUE, xlab="log2(intensity)",
                          trFunc=log2xplus1)}
\arguments{
   \item{MSnSetObj}{MSnSet; an object of class MSnSet}
   \item{sampleColours}{character: a vector of colors for samples}
   \item{title}{character: title for the plot}
   \item{colourBy}{character: column name from pData(MSnSetObj) to use for coloring samples}
   \item{transform}{logical: whether to log transform intensities}
   \item{xlab}{character: label for x-axis}
   \item{trFunc}{func: internal helper function for log transformation}
}
\details{
  The column provided to the "colourBy" argument will be used to colour the samples. The colours
  will be determined using the function \link{assignColours}, alternatively the user may specify a named
  vector of colours using the "sampleColours" argument. The names of the "sampleColours" vector should
  match the unique values in the "colourBy" column.
}
\value{
An intensity distribution plot for quantitative proteomics data.
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER")

# custom colours
customCols <- rainbow(length(unique(pData(MSnSet_data)$SampleGroup)))
names(customCols) <- unique(pData(MSnSet_data)$SampleGroup)
intensityPlot(MSnSet_data, title = "qPLEX_RIME_ER", sampleColours = customCols)
}
