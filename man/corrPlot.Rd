\name{corrPlot}
\alias{corrPlot}
\title{Correlation plot}
\description{
  Computes and display correlation plot for samples within MSnSet
}
\usage{corrPlot(MSnSetObj, addValues=TRUE, title="")}
\arguments{
   \item{MSnSetObj}{MSnSet; an object of class MSnSet}
   \item{addValues}{logical: adds correlation values to the plot}
   \item{title}{character; title of the plot}  
}
\value{(Invisibly) returns a reordered correlation matrix plot.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),Sequences=2,Accessions=6)
corrPlot(MSnSet_data,addValues=TRUE,title="Correlation plot")
}
