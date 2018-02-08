\name{hierarchicalPlot}
\alias{hierarchicalPlot}
\title{Hierarchical plot}
\description{
  Computes and display hierarchical clustering plot for samples within MSnSet
}
\usage{hierarchicalPlot(data,label_color,branchlength=20,title)}
\arguments{
   \item{data}{MSnSet; An object of class MSnSet}
   \item{label_color}{Character: A vector of colors for samples}
   \item{branchlength}{numeric: How long to make the colored part of the branches}
   \item{title}{Character: The main title for the dendrogram}
}
\examples{
data(human_anno)
data(exp3_OHT_ESR1)
MSnset_data <- convertToMSnset(exp3_OHT_ESR1$intensities_qPLEX1,metadata=exp3_OHT_ESR1$metadata_qPLEX1,
indExpData=c(7:16),indFData=c(1:6))
exprs(MSnset_data) <- exprs(MSnset_data)+0.01
label_color <- c(rep(c("red","darkgreen","yellow","brown"),2),rep("orange",2))
hierarchicalPlot(MSnset_data,label_color,branchlength=55,title="qPLEX_RIME_ER")
}
