\name{plotMeanVar}
\alias{plotMeanVar}
\title{Mean variance plot}
\description{
  Computes and plots mean-variance for samples in MSnset
}
\usage{plotMeanVar(MSnSetObj, title="")}
\arguments{
   \item{MSnSetObj}{MSnSet; an object of class MSnSet}
   \item{title}{character: title for the plot}
}
\value{Plot of mean variance.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnSet_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1, metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16), Sequences=2, Accessions=6)
plotMeanVar(MSnSet_data, title="Mean_Variance")
}
