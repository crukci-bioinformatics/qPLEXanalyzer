\name{intensityPlot}
\alias{intensityPlot}
\title{Intensity Distribution Plot}
\description{
  Intensity distribution plot of all the samples}
\usage{intensityPlot(data, sampleColours, log2Transform=TRUE, title = "", 
                                      xlab = "log2(intensity)", minIntensity = NA, 
                                      maxIntensity = NA)}
\arguments{
   \item{data}{MSnSet; An object of class MSnSet}
   \item{sampleColours}{Character: A named vector that maps samples to colours
                                    generated from function assignColours}
   \item{log2Transform}{Logical: Whether to log transform intensities}
   \item{title}{Character: Title for the plot}
   \item{xlab}{Character: Text for x-axis}
   \item{minIntensity}{numeric: minimum intensity to be plotted}
   \item{maxIntensity}{numeric: maximum intensity to be plotted}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
sampleColours <- assignColours(MSnset_data, groupColourSchemes = c("Blues","Greens","Purples",
"Reds","YlOrBr"))
print(intensityPlot(MSnset_data, sampleColours, title = "qPLEX_RIME_ER"))
}
