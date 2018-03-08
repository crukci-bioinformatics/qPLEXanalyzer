\name{hierarchicalPlot}
\alias{hierarchicalPlot}
\title{Hierarchical plot}
\description{
  Computes and display hierarchical clustering plot for samples within MSnSet
}
\usage{hierarchicalPlot(MSnSetObj,sampleColours=NULL,colourBy="SampleGroup", horizontal=TRUE,
                             title="")}
\arguments{
   \item{MSnSetObj}{MSnSet; An object of class MSnSet}
   \item{sampleColours}{Character: A vector of colors for samples}
   \item{colourBy}{Character: column name to use for coloring samples from pData(MSnSetObj)}
   \item{horizontal}{Logical: Define orientation of the dendrogram"}
   \item{title}{Character: The main title for the dendrogram}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
exprs(MSnSet_data) <- exprs(MSnSet_data)+0.01
hierarchicalPlot(MSnSet_data, title="qPLEX_RIME_ER")
}
