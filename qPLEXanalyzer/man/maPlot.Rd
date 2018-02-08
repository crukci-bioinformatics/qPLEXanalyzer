\name{maPlot}
\alias{maPlot}
\title{MA Plot}
\description{
  MA plot of differential statistics results
}
\usage{maPlot(differentialExpressionResults, selectedGenes = NULL, xlab = "average log2(intensity)",
                   ylab = "log2FC", significanceLevel = 0.05, title="",
                   minLogFoldChangeForLabelling = 1.5, controlLogFoldChangeThreshold = -Inf, 
                   pointSize = 1.5)}
\arguments{
   \item{differentialExpressionResults}{data frame; Output of the function getContrastResults}
   \item{selectedGenes}{Character: A vector of gene names to label differently}
   \item{xlab}{Character: text for x-axis}
   \item{ylab}{Character: text for y-axis}
   \item{significanceLevel}{numeric: threshold for p-value display}
   \item{title}{Character: title for the plot}
   \item{minLogFoldChangeForLabelling}{numeric: minimum threshold for labeling the fold change}
   \item{controlLogFoldChangeThreshold}{numeric: plot genes above controlLogFoldChange threshold}
   \item{pointSize}{numeric: size of points on plot}
}
\examples{
data(human_anno)
data(results)
maPlot(results, significanceLevel = 0.05, minLogFoldChangeForLabelling = 1.0)
}
