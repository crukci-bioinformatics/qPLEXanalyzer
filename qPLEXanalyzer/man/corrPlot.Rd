\name{corrPlot}
\alias{corrPlot}
\title{Correlation plot}
\description{
  Computes and display correlation plot for samples within MSnSet
}
\usage{corrPlot(data,method="shade",title="Correlation plot")}
\arguments{
   \item{data}{MSnSet; An object of class MSnSet}
   \item{method}{Character: the visualization method for displaying correlation plot. 
      It supports seven methods, "circle", "square", "ellipse", "number",
      "pie", "shade" and "color"}
   \item{title}{character; title of the plot}  
}
\value{(Invisibly) returns a reordered correlation matrix.}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
corrPlot(MSnset_data,method="shade")
}
