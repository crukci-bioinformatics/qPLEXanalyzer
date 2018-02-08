\name{plotMeanVar}
\alias{plotMeanVar}
\title{Mean variance plot}
\description{
  Computes and plots mean-variance for samples within MSnset
}
\usage{plotMeanVar(data,title,color,ylim=NULL,xlim=NULL)}
\arguments{
   \item{data}{MSnSet; An object of class MSnSet}
   \item{title}{Character: title for the plot}
   \item{color}{Character: colors for plotting}
   \item{ylim}{numeric}
   \item{xlim}{numeric}
}
\value{Plot of mean variance.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,indExpData=c(7:16),indFData=c(1:6))
plotMeanVar(MSnset_data,title="Mean_Variance",color="red")
}
