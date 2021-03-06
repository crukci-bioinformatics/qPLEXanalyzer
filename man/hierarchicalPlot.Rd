% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hierarchicalPlot.R
\name{hierarchicalPlot}
\alias{hierarchicalPlot}
\title{Hierarchical clustering plot}
\usage{
hierarchicalPlot(
  MSnSetObj,
  sampleColours = NULL,
  colourBy = "SampleGroup",
  horizontal = TRUE,
  title = ""
)
}
\arguments{
\item{MSnSetObj}{MSnSet; an object of class MSnSet}

\item{sampleColours}{character: a named vector of colors for samples, names
should be values of \code{colourBy} column}

\item{colourBy}{character: column name from \code{pData(MSnSetObj)} to use 
for coloring samples}

\item{horizontal}{logical: define orientation of the dendrogram}

\item{title}{character: the main title for the dendrogram}
}
\value{
An object created by \code{ggplot}
}
\description{
Computes and displays hierarchical clustering plot for samples in MSnSet
}
\examples{

data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, 
                               metadata=exp3_OHT_ESR1$metadata_qPLEX1,
                               indExpData=c(7:16), 
                               Sequences=2, 
                               Accessions=6)
exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
hierarchicalPlot(MSnSet_data, title="qPLEX_RIME_ER")

}
