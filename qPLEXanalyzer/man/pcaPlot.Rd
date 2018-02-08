\name{pcaPlot}
\alias{pcaPlot}
\title{PCA plot}
\description{
  PCA plots of the samples within MSnset
}
\usage{pcaPlot(data, groupColours, title = "", labels = NULL, legend = TRUE, labelsize=2.5, logTransform = TRUE)}
\arguments{
   \item{data}{MSnSet; An object of class MSnSet}
  \item{groupColours}{Character: vector of colors for samples}
  \item{title}{Character: title for the plot}
  \item{labels}{Character: A vector of labels for samples}
  \item{legend}{Logical: Whether to display legends or not}
  \item{labelsize}{Numeric: Size of the labels}
  \item{logTransform}{Logical: Whether to logTransform intensities}
}
\value{An instance of class "MSnSet".}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
exprs(MSnset_data) <- exprs(MSnset_data)+0.01
grpcolor <- c("red","darkgreen","yellow","brown","orange")
pcaPlot(MSnset_data, groupColours=grpcolor, title = "MSnset_data", labels = MSnset_data$Experiment,
labelsize=2.5, legend=TRUE)
}
